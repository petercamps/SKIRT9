/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "H5Core.hpp"
#include "FatalError.hpp"

#ifdef BUILD_HDF5
#    include "StringUtils.hpp"
#    include "System.hpp"
#    include <hdf5.h>
#    include <chrono>
#    include <condition_variable>
#    include <mutex>
#    include <set>
#endif

////////////////////////////////////////////////////////////////////

#ifdef BUILD_HDF5

namespace
{
    // the handles in the interface of this namespace must be able to hold the identifiers used by the HDF5 library
    static_assert(sizeof(hid_t) <= sizeof(std::int64_t), "HDF5 identifiers do not fit in a 64-bit integer");

    // the mutex that protects all calls into the HDF5 library
    std::mutex _mutex;

    // the condition that is notified when an object is released, for functions waiting until a file is no longer
    // held open for reading
    std::condition_variable _released;

    // the time a function waits until a file is no longer held open for reading before it reports an error
    constexpr auto _readerTimeLimit = std::chrono::minutes(10);

    // An instance of this class locks the mutex for its lifetime. It also stops the HDF5 library from printing its
    // own error messages for the calling thread, because failures are reported through fatal errors instead.
    class Lock
    {
        std::unique_lock<std::mutex> _lock;

    public:
        Lock() : _lock(_mutex) { H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr); }

        // releases the mutex until an object is released or the specified time has passed, and then locks it again
        void wait(std::chrono::milliseconds duration) { _released.wait_for(_lock, duration); }
    };

    // An instance of this class holds an identifier handed out by the HDF5 library, and releases it upon
    // destruction through the function specified at construction, unless the identifier is invalid.
    class Id
    {
        hid_t _id{-1};
        herr_t (*_close)(hid_t){nullptr};

    public:
        Id(hid_t id, herr_t (*close)(hid_t)) : _id(id), _close(close) {}
        ~Id()
        {
            if (_id >= 0) _close(_id);
        }
        Id(const Id&) = delete;
        Id& operator=(const Id&) = delete;

        // returns the identifier
        hid_t get() const { return _id; }

        // returns true if the identifier is valid
        bool valid() const { return _id >= 0; }

        // returns the identifier and gives up the responsibility for releasing it
        hid_t release()
        {
            hid_t id = _id;
            _id = -1;
            return id;
        }
    };

    // opens the HDF5 file with the specified path for reading and returns its identifier, which must be released by
    // the caller; throws a fatal error if the file does not exist or cannot be opened as an HDF5 file
    hid_t openFile(const string& filePath)
    {
        // on a file system that does not support file locking, proceed without locking
        Id access(H5Pcreate(H5P_FILE_ACCESS), H5Pclose);
#    if H5_VERSION_GE(1, 12, 1) || (H5_VERSION_GE(1, 10, 7) && !H5_VERSION_GE(1, 12, 0))
        H5Pset_file_locking(access.get(), true, true);
#    endif
        Id file(H5Fopen(filePath.c_str(), H5F_ACC_RDONLY, access.get()), H5Fclose);
        if (!file.valid())
        {
            if (!System::isFile(filePath)) throw FATALERROR("HDF5 file does not exist: " + filePath);
            throw FATALERROR("File cannot be opened as an HDF5 file: " + filePath);
        }
        return file.release();
    }

    // opens the attribute with the specified name of the specified object, verifies that it holds a single value,
    // and returns its identifier, which must be released by the caller
    hid_t openScalarAttribute(hid_t object, const string& name, const string& context)
    {
        Id attribute(H5Aopen(object, name.c_str(), H5P_DEFAULT), H5Aclose);
        if (!attribute.valid()) throw FATALERROR("Missing attribute '" + name + "' for " + context);

        Id space(H5Aget_space(attribute.get()), H5Sclose);
        if (H5Sget_simple_extent_npoints(space.get()) != 1)
            throw FATALERROR("Attribute '" + name + "' for " + context + " does not hold a single value");
        return attribute.release();
    }

    // returns the class of the data type of the specified attribute
    H5T_class_t attributeClass(hid_t attribute)
    {
        Id type(H5Aget_type(attribute), H5Tclose);
        return H5Tget_class(type.get());
    }

    // returns the class of the data type of the specified dataset
    H5T_class_t datasetClass(hid_t dataset)
    {
        Id type(H5Dget_type(dataset), H5Tclose);
        return H5Tget_class(type.get());
    }

    // returns the number of elements in the specified dataset
    size_t numElements(hid_t dataset, const string& context)
    {
        Id space(H5Dget_space(dataset), H5Sclose);
        hssize_t count = H5Sget_simple_extent_npoints(space.get());
        if (count < 0) throw FATALERROR("Cannot determine the size of " + context);
        return static_cast<size_t>(count);
    }

    // This function is called by the HDF5 library when a value cannot be converted without loss of information. It
    // sets the Boolean flag passed as user data, and causes the read operation to fail.
    H5T_conv_ret_t reportConversionFailure(H5T_conv_except_t /*exceptType*/, hid_t /*sourceType*/, hid_t /*destType*/,
                                           void* /*sourceBuffer*/, void* /*destBuffer*/, void* userData)
    {
        *static_cast<bool*>(userData) = true;
        return H5T_CONV_ABORT;
    }

    // reads the complete contents of the specified dataset into the specified buffer, converting each value to the
    // specified memory data type; if requested, the function throws a fatal error when a value cannot be converted
    // without loss of information, and otherwise such a value is clipped to the range of the memory data type
    void readValues(hid_t dataset, hid_t memoryType, void* buffer, bool lossless, const string& context)
    {
        bool lost = false;
        Id transfer(H5Pcreate(H5P_DATASET_XFER), H5Pclose);
        if (lossless) H5Pset_type_conv_cb(transfer.get(), reportConversionFailure, &lost);

        herr_t status = H5Dread(dataset, memoryType, H5S_ALL, H5S_ALL, transfer.get(), buffer);
        if (lost) throw FATALERROR("Some values in " + context + " cannot be represented without loss of information");
        if (status < 0) throw FATALERROR("Cannot read the values in " + context);
    }
}

////////////////////////////////////////////////////////////////////

bool H5Core::available()
{
    return true;
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::openGroup(const string& filePath, const string& groupPath)
{
    Lock lock;

    // open the file and the group; the library keeps the file open as long as the group is open
    Id file(openFile(filePath), H5Fclose);
    Id group(H5Gopen2(file.get(), ("/" + groupPath).c_str(), H5P_DEFAULT), H5Gclose);
    if (!group.valid()) throw FATALERROR("HDF5 file " + filePath + " has no group '" + groupPath + "'");
    return group.release();
}

////////////////////////////////////////////////////////////////////

bool H5Core::hasGroup(const string& filePath, const string& groupPath)
{
    Lock lock;

    if (!System::isFile(filePath)) return false;
    Id file(openFile(filePath), H5Fclose);
    Id group(H5Gopen2(file.get(), ("/" + groupPath).c_str(), H5P_DEFAULT), H5Gclose);
    return group.valid();
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::openGroup(std::int64_t group, const string& name, const string& context)
{
    Lock lock;

    Id child(H5Gopen2(group, name.c_str(), H5P_DEFAULT), H5Gclose);
    if (!child.valid()) throw FATALERROR("There is no group '" + name + "' in " + context);
    return child.release();
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the names of the objects of the specified type that are located directly inside the specified group,
    // in alphabetical order; the caller must hold the lock
    vector<string> memberNames(hid_t group, H5I_type_t type, const string& context)
    {
        H5G_info_t info;
        if (H5Gget_info(group, &info) < 0) throw FATALERROR("Cannot list the contents of " + context);

        vector<string> names;
        for (hsize_t index = 0; index != info.nlinks; ++index)
        {
            // get the name of the link with this index, in two steps because the length of the name is not known
            ssize_t length = H5Lget_name_by_idx(group, ".", H5_INDEX_NAME, H5_ITER_INC, index, nullptr, 0, H5P_DEFAULT);
            if (length < 0) throw FATALERROR("Cannot list the contents of " + context);
            string name(static_cast<size_t>(length) + 1, '\0');  // includes room for the terminating null character
            H5Lget_name_by_idx(group, ".", H5_INDEX_NAME, H5_ITER_INC, index, &name[0], name.size(), H5P_DEFAULT);
            name.resize(static_cast<size_t>(length));

            // include the name only if it refers to an object of the requested type
            Id object(H5Oopen(group, name.c_str(), H5P_DEFAULT), H5Oclose);
            if (object.valid() && H5Iget_type(object.get()) == type) names.push_back(name);
        }
        return names;
    }
}

////////////////////////////////////////////////////////////////////

vector<string> H5Core::groupNames(std::int64_t group, const string& context)
{
    Lock lock;
    return memberNames(group, H5I_GROUP, context);
}

////////////////////////////////////////////////////////////////////

vector<string> H5Core::datasetNames(std::int64_t group, const string& context)
{
    Lock lock;
    return memberNames(group, H5I_DATASET, context);
}

////////////////////////////////////////////////////////////////////

bool H5Core::hasDataset(std::int64_t group, const string& name)
{
    Lock lock;

    if (H5Lexists(group, name.c_str(), H5P_DEFAULT) <= 0) return false;
    Id object(H5Oopen(group, name.c_str(), H5P_DEFAULT), H5Oclose);
    return object.valid() && H5Iget_type(object.get()) == H5I_DATASET;
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::openDataset(std::int64_t group, const string& name, const string& context)
{
    Lock lock;

    Id dataset(H5Dopen2(group, name.c_str(), H5P_DEFAULT), H5Dclose);
    if (!dataset.valid()) throw FATALERROR("There is no dataset '" + name + "' in " + context);
    return dataset.release();
}

////////////////////////////////////////////////////////////////////

void H5Core::close(std::int64_t object) noexcept
{
    if (object >= 0)
    {
        try
        {
            Lock lock;
            H5Oclose(object);
        }
        catch (...)
        {}
        _released.notify_all();
    }
}

////////////////////////////////////////////////////////////////////

bool H5Core::hasAttribute(std::int64_t object, const string& name)
{
    Lock lock;

    return H5Aexists(object, name.c_str()) > 0;
}

////////////////////////////////////////////////////////////////////

bool H5Core::boolAttribute(std::int64_t object, const string& name, const string& context)
{
    Lock lock;

    Id attribute(openScalarAttribute(object, name, context), H5Aclose);
    H5T_class_t type = attributeClass(attribute.get());
    if (type != H5T_INTEGER && type != H5T_ENUM)
        throw FATALERROR("Attribute '" + name + "' for " + context + " is not a Boolean value");

    // a value outside of the range is clipped by the library, so that a nonzero value remains nonzero
    long long value = 0;
    if (H5Aread(attribute.get(), H5T_NATIVE_LLONG, &value) < 0)
        throw FATALERROR("Cannot read attribute '" + name + "' for " + context);
    return value != 0;
}

////////////////////////////////////////////////////////////////////

int H5Core::intAttribute(std::int64_t object, const string& name, const string& context)
{
    Lock lock;

    Id attribute(openScalarAttribute(object, name, context), H5Aclose);
    H5T_class_t type = attributeClass(attribute.get());

    // an integer value must be in range; a value outside of the 64-bit range is clipped by the library
    if (type == H5T_INTEGER)
    {
        long long value = 0;
        if (H5Aread(attribute.get(), H5T_NATIVE_LLONG, &value) < 0)
            throw FATALERROR("Cannot read attribute '" + name + "' for " + context);
        if (value < INT_MIN || value > INT_MAX)
            throw FATALERROR("Attribute '" + name + "' for " + context + " does not fit in a 32-bit integer");
        return static_cast<int>(value);
    }

    // a floating point value must be in range and have no fractional part
    if (type == H5T_FLOAT)
    {
        double value = 0.;
        if (H5Aread(attribute.get(), H5T_NATIVE_DOUBLE, &value) < 0)
            throw FATALERROR("Cannot read attribute '" + name + "' for " + context);
        if (!(value >= INT_MIN && value <= INT_MAX) || value != std::trunc(value))
            throw FATALERROR("Attribute '" + name + "' for " + context + " is not a 32-bit integer");
        return static_cast<int>(value);
    }

    throw FATALERROR("Attribute '" + name + "' for " + context + " is not an integer");
}

////////////////////////////////////////////////////////////////////

double H5Core::doubleAttribute(std::int64_t object, const string& name, const string& context)
{
    Lock lock;

    Id attribute(openScalarAttribute(object, name, context), H5Aclose);
    H5T_class_t type = attributeClass(attribute.get());
    if (type != H5T_INTEGER && type != H5T_FLOAT)
        throw FATALERROR("Attribute '" + name + "' for " + context + " is not a number");

    double value = 0.;
    if (H5Aread(attribute.get(), H5T_NATIVE_DOUBLE, &value) < 0)
        throw FATALERROR("Cannot read attribute '" + name + "' for " + context);
    return value;
}

////////////////////////////////////////////////////////////////////

string H5Core::stringAttribute(std::int64_t object, const string& name, const string& context)
{
    Lock lock;

    Id attribute(openScalarAttribute(object, name, context), H5Aclose);
    Id type(H5Aget_type(attribute.get()), H5Tclose);
    if (H5Tget_class(type.get()) != H5T_STRING)
        throw FATALERROR("Attribute '" + name + "' for " + context + " is not a string");

    // the memory data type is a C string with the character set of the stored string (ASCII or UTF-8)
    Id memoryType(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_cset(memoryType.get(), H5Tget_cset(type.get()));

    // a string of variable length is allocated by the library, and must be released after copying it
    if (H5Tis_variable_str(type.get()) > 0)
    {
        H5Tset_size(memoryType.get(), H5T_VARIABLE);
        char* buffer = nullptr;
        if (H5Aread(attribute.get(), memoryType.get(), &buffer) < 0)
            throw FATALERROR("Cannot read attribute '" + name + "' for " + context);
        string value = buffer ? buffer : "";
        H5free_memory(buffer);
        return value;
    }

    // a string of fixed length is converted to a null-terminated string, which removes any padding
    size_t size = H5Tget_size(type.get());
    H5Tset_size(memoryType.get(), size + 1);
    H5Tset_strpad(memoryType.get(), H5T_STR_NULLTERM);
    string value(size + 1, '\0');
    if (H5Aread(attribute.get(), memoryType.get(), &value[0]) < 0)
        throw FATALERROR("Cannot read attribute '" + name + "' for " + context);
    value.resize(value.find('\0'));
    return value;
}

////////////////////////////////////////////////////////////////////

vector<size_t> H5Core::shape(std::int64_t dataset, const string& context)
{
    Lock lock;

    Id space(H5Dget_space(dataset), H5Sclose);
    int rank = H5Sget_simple_extent_ndims(space.get());
    if (rank < 0) throw FATALERROR("Cannot determine the shape of " + context);

    // a dataset without any data is reported as an empty dataset with a single dimension
    if (H5Sget_simple_extent_type(space.get()) == H5S_NULL) return vector<size_t>{0};

    vector<hsize_t> dims(static_cast<size_t>(rank));
    if (rank > 0) H5Sget_simple_extent_dims(space.get(), dims.data(), nullptr);
    return vector<size_t>(dims.begin(), dims.end());
}

////////////////////////////////////////////////////////////////////

vector<char> H5Core::readBool(std::int64_t dataset, const string& context)
{
    Lock lock;

    H5T_class_t type = datasetClass(dataset);
    if (type != H5T_INTEGER && type != H5T_ENUM) throw FATALERROR("The values in " + context + " are not Booleans");

    // read the values as signed bytes: a value outside of that range is clipped by the library, so that a nonzero
    // value remains nonzero (also if it is negative), and there is no need to require lossless conversion; the
    // library stores the bytes directly into the result, where a signed byte occupies the same storage as a char
    vector<char> values(numElements(dataset, context));
    if (!values.empty()) readValues(dataset, H5T_NATIVE_SCHAR, values.data(), false, context);

    // reduce each value to one for true or zero for false
    for (char& value : values) value = value != 0;
    return values;
}

////////////////////////////////////////////////////////////////////

vector<int> H5Core::readInt(std::int64_t dataset, const string& context)
{
    Lock lock;

    H5T_class_t type = datasetClass(dataset);
    if (type != H5T_INTEGER && type != H5T_FLOAT) throw FATALERROR("The values in " + context + " are not integers");

    vector<int> values(numElements(dataset, context));
    if (!values.empty()) readValues(dataset, H5T_NATIVE_INT, values.data(), true, context);
    return values;
}

////////////////////////////////////////////////////////////////////

Array H5Core::readDouble(std::int64_t dataset, const string& context)
{
    Lock lock;

    H5T_class_t type = datasetClass(dataset);
    if (type != H5T_INTEGER && type != H5T_FLOAT) throw FATALERROR("The values in " + context + " are not numbers");

    Array values(numElements(dataset, context));
    if (values.size()) readValues(dataset, H5T_NATIVE_DOUBLE, &values[0], true, context);
    return values;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // opens the HDF5 file with the specified path for writing, creating it if it does not exist, and returns its
    // identifier, which must be released by the caller; the caller must hold the lock
    // returns the name of the specified open file, as passed to the function that opened it
    string fileName(hid_t file)
    {
        ssize_t length = H5Fget_name(file, nullptr, 0);
        if (length <= 0) return string();
        string name(length + 1, '\0');
        H5Fget_name(file, &name[0], name.size());
        name.resize(length);
        return name;
    }

    // calls the specified function with the name and the access flags of the file of each object held open by this
    // process; the caller must hold the lock
    template<class Function> void forEachOpenObject(Function function)
    {
        ssize_t count = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
        if (count <= 0) return;
        vector<hid_t> objects(count);
        count = H5Fget_obj_ids(H5F_OBJ_ALL, H5F_OBJ_ALL, objects.size(), objects.data());
        for (ssize_t index = 0; index < count; ++index)
        {
            Id file(H5Iget_file_id(objects[index]), H5Fclose);
            unsigned int intent = 0;
            if (file.valid() && H5Fget_intent(file.get(), &intent) >= 0) function(fileName(file.get()), intent);
        }
    }

    // returns true if an object held open by this process belongs to the specified file, opened for reading only;
    // the caller must hold the lock
    bool isHeldForReading(const string& filePath)
    {
        if (!System::isFile(filePath)) return false;
        string target = System::canonicalPath(filePath);
        bool held = false;
        forEachOpenObject([&held, &target](const string& name, unsigned int intent) {
            if (!(intent & H5F_ACC_RDWR) && System::canonicalPath(name) == target) held = true;
        });
        return held;
    }

    // waits until no object held open by this process belongs to the specified file opened for reading only, because
    // the library cannot reopen such a file for writing; throws a fatal error if this takes longer than the time limit
    void waitForReaders(Lock& lock, const string& filePath)
    {
        auto deadline = std::chrono::steady_clock::now() + _readerTimeLimit;
        while (isHeldForReading(filePath))
        {
            if (std::chrono::steady_clock::now() > deadline)
                throw FATALERROR("HDF5 file remains open for reading in this process, so that it cannot be opened "
                                 "for writing: "
                                 + filePath);
            lock.wait(std::chrono::milliseconds(100));
        }
    }

    hid_t openFileForWriting(const string& filePath)
    {
        // on a file system that does not support file locking, proceed without locking
        Id access(H5Pcreate(H5P_FILE_ACCESS), H5Pclose);
#    if H5_VERSION_GE(1, 12, 1) || (H5_VERSION_GE(1, 10, 7) && !H5_VERSION_GE(1, 12, 0))
        H5Pset_file_locking(access.get(), true, true);
#    endif

        // open an existing file
        if (System::isFile(filePath))
        {
            Id file(H5Fopen(filePath.c_str(), H5F_ACC_RDWR, access.get()), H5Fclose);
            if (!file.valid())
            {
#    if H5_VERSION_GE(1, 12, 0)
                bool isHdf5 = H5Fis_accessible(filePath.c_str(), access.get()) > 0;
#    else
                bool isHdf5 = H5Fis_hdf5(filePath.c_str()) > 0;
#    endif
                if (!isHdf5) throw FATALERROR("File cannot be opened as an HDF5 file: " + filePath);
                throw FATALERROR("HDF5 file cannot be opened for writing: " + filePath
                                 + "\n(it may be read-only, still open for reading, or in use by another process)");
            }
            return file.release();
        }

        // or create a new file
        Id file(H5Fcreate(filePath.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, access.get()), H5Fclose);
        if (!file.valid()) throw FATALERROR("Cannot create HDF5 file: " + filePath);
        return file.release();
    }

    // returns the type of the object at the specified path relative to the specified location, or H5I_BADID if there
    // is no such object; the caller must hold the lock
    H5I_type_t objectType(hid_t location, const string& path)
    {
        // verify that each segment of the path exists, because the library does not allow testing a nested path at once
        string partial;
        for (const string& segment : StringUtils::split(path, "/"))
        {
            partial += (partial.empty() ? "" : "/") + segment;
            if (H5Lexists(location, partial.c_str(), H5P_DEFAULT) <= 0) return H5I_BADID;
        }
        Id object(H5Oopen(location, path.c_str(), H5P_DEFAULT), H5Oclose);
        return object.valid() ? H5Iget_type(object.get()) : H5I_BADID;
    }

    // returns true if the specified group contains at least one object of the specified type; if the last argument is
    // true, only the first member of the group (in the order of the names) is inspected, which suffices as a sanity
    // check and remains fast for a group with many members; the caller must hold the lock
    bool containsObjectOfType(hid_t group, H5I_type_t type, bool firstOnly, const string& context)
    {
        H5G_info_t info;
        if (H5Gget_info(group, &info) < 0) throw FATALERROR("Cannot list the contents of " + context);
        hsize_t count = firstOnly ? std::min(info.nlinks, static_cast<hsize_t>(1)) : info.nlinks;
        for (hsize_t index = 0; index != count; ++index)
        {
            Id object(H5Oopen_by_idx(group, ".", H5_INDEX_NAME, H5_ITER_INC, index, H5P_DEFAULT), H5Oclose);
            if (object.valid() && H5Iget_type(object.get()) == type) return true;
        }
        return false;
    }

    // returns the identifier of a new data type for Booleans, in the representation used by the h5py package: an
    // enumeration based on a signed byte with the members FALSE (0) and TRUE (1); the caller must hold the lock
    hid_t booleanType()
    {
        hid_t type = H5Tenum_create(H5T_NATIVE_SCHAR);
        signed char value = 0;
        H5Tenum_insert(type, "FALSE", &value);
        value = 1;
        H5Tenum_insert(type, "TRUE", &value);
        return type;
    }

    // returns the identifier of a new data type for strings of variable length in UTF-8 encoding; the caller must
    // hold the lock
    hid_t stringType()
    {
        hid_t type = H5Tcopy(H5T_C_S1);
        H5Tset_size(type, H5T_VARIABLE);
        H5Tset_cset(type, H5T_CSET_UTF8);
        return type;
    }

    // writes an attribute with the specified name, file and memory data types, and value to the specified object,
    // replacing an existing attribute with the same name; the caller must hold the lock
    void writeAttribute(hid_t object, const string& name, hid_t fileType, hid_t memoryType, const void* value,
                        const string& context)
    {
        if (H5Aexists(object, name.c_str()) > 0 && H5Adelete(object, name.c_str()) < 0)
            throw FATALERROR("Cannot replace attribute '" + name + "' for " + context);
        Id space(H5Screate(H5S_SCALAR), H5Sclose);
        Id attribute(H5Acreate2(object, name.c_str(), fileType, space.get(), H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
        if (!attribute.valid() || H5Awrite(attribute.get(), memoryType, value) < 0)
            throw FATALERROR("Cannot write attribute '" + name + "' for " + context);
    }

    // creates a dataset with the specified name, shape, and file data type in the specified group, writes the
    // specified values in the specified memory data type to it, and returns its identifier, which must be released by
    // the caller; the caller must hold the lock
    hid_t writeDataset(hid_t group, const string& name, const vector<size_t>& shape, hid_t fileType, hid_t memoryType,
                       const void* values, const string& context)
    {
        vector<hsize_t> dims(shape.begin(), shape.end());
        Id space(H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr), H5Sclose);
        if (!space.valid()) throw FATALERROR("Cannot create dataset '" + name + "' in " + context);
        Id dataset(H5Dcreate2(group, name.c_str(), fileType, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                   H5Dclose);
        if (!dataset.valid()) throw FATALERROR("Cannot create dataset '" + name + "' in " + context);
        if (H5Sget_simple_extent_npoints(space.get()) > 0
            && H5Dwrite(dataset.get(), memoryType, H5S_ALL, H5S_ALL, H5P_DEFAULT, values) < 0)
            throw FATALERROR("Cannot write dataset '" + name + "' in " + context);
        return dataset.release();
    }
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::createGroup(const string& filePath, const string& groupPath)
{
    Lock lock;
    waitForReaders(lock, filePath);

    Id file(openFileForWriting(filePath), H5Fclose);
    string location = filePath + ":" + groupPath;

    // each parent of the new group must be a suite: a group that does not directly contain datasets
    auto segments = StringUtils::split(groupPath, "/");
    string parentPath;
    for (size_t index = 0; index + 1 < segments.size(); ++index)
    {
        parentPath += (parentPath.empty() ? "" : "/") + segments[index];
        H5I_type_t type = objectType(file.get(), parentPath);
        if (type == H5I_BADID) break;
        bool isSuite = type == H5I_GROUP;
        if (isSuite)
        {
            Id parent(H5Gopen2(file.get(), parentPath.c_str(), H5P_DEFAULT), H5Gclose);
            isSuite = !containsObjectOfType(parent.get(), H5I_DATASET, true, "HDF5 group " + parentPath);
        }
        if (!isSuite)
            throw FATALERROR("Cannot create bundle " + location + " because '" + parentPath + "' is not a suite");
    }

    // an existing object with the same path must be a bundle: a group that does not contain other groups
    H5I_type_t type = objectType(file.get(), groupPath);
    if (type != H5I_BADID)
    {
        if (type != H5I_GROUP)
            throw FATALERROR("Cannot create bundle " + location + " because there is a dataset with that path");
        {
            Id existing(H5Gopen2(file.get(), groupPath.c_str(), H5P_DEFAULT), H5Gclose);
            if (containsObjectOfType(existing.get(), H5I_GROUP, false, "HDF5 group " + location))
                throw FATALERROR("Cannot create bundle " + location + " because there is a suite with that path");
        }
        if (H5Ldelete(file.get(), groupPath.c_str(), H5P_DEFAULT) < 0)
            throw FATALERROR("Cannot replace the existing bundle " + location);
    }

    // create the group, and any missing parent groups
    Id links(H5Pcreate(H5P_LINK_CREATE), H5Pclose);
    H5Pset_create_intermediate_group(links.get(), 1);
    Id group(H5Gcreate2(file.get(), groupPath.c_str(), links.get(), H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
    if (!group.valid()) throw FATALERROR("Cannot create bundle " + location);
    return group.release();
}

////////////////////////////////////////////////////////////////////

void H5Core::setBoolAttribute(std::int64_t object, const string& name, bool value, const string& context)
{
    Lock lock;

    Id type(booleanType(), H5Tclose);
    signed char flag = value ? 1 : 0;
    writeAttribute(object, name, type.get(), type.get(), &flag, context);
}

////////////////////////////////////////////////////////////////////

void H5Core::setIntAttribute(std::int64_t object, const string& name, int value, const string& context)
{
    Lock lock;

    writeAttribute(object, name, H5T_STD_I32LE, H5T_NATIVE_INT, &value, context);
}

////////////////////////////////////////////////////////////////////

void H5Core::setDoubleAttribute(std::int64_t object, const string& name, double value, const string& context)
{
    Lock lock;

    writeAttribute(object, name, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, &value, context);
}

////////////////////////////////////////////////////////////////////

void H5Core::setStringAttribute(std::int64_t object, const string& name, const string& value, const string& context)
{
    Lock lock;

    Id type(stringType(), H5Tclose);
    const char* buffer = value.c_str();
    writeAttribute(object, name, type.get(), type.get(), &buffer, context);
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::createBoolDataset(std::int64_t group, const string& name, const vector<size_t>& shape,
                                       const char* values, const string& context)
{
    Lock lock;

    // a character occupies the same storage as a signed byte, the base type of the enumeration
    Id type(booleanType(), H5Tclose);
    return writeDataset(group, name, shape, type.get(), type.get(), values, context);
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::createIntDataset(std::int64_t group, const string& name, const vector<size_t>& shape,
                                      const int* values, const string& context)
{
    Lock lock;

    return writeDataset(group, name, shape, H5T_STD_I32LE, H5T_NATIVE_INT, values, context);
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::createDoubleDataset(std::int64_t group, const string& name, const vector<size_t>& shape,
                                         const double* values, const string& context)
{
    Lock lock;

    return writeDataset(group, name, shape, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, values, context);
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::createTextDataset(std::int64_t group, const string& name, const string& text,
                                       const string& context)
{
    Lock lock;

    Id type(stringType(), H5Tclose);
    Id space(H5Screate(H5S_SCALAR), H5Sclose);
    Id dataset(H5Dcreate2(group, name.c_str(), type.get(), space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
               H5Dclose);
    if (!dataset.valid()) throw FATALERROR("Cannot create dataset '" + name + "' in " + context);
    const char* buffer = text.c_str();
    if (H5Dwrite(dataset.get(), type.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer) < 0)
        throw FATALERROR("Cannot write dataset '" + name + "' in " + context);
    return dataset.release();
}

////////////////////////////////////////////////////////////////////

void H5Core::verifyWritable(const string& filePath)
{
    Lock lock;
    waitForReaders(lock, filePath);
    Id file(openFileForWriting(filePath), H5Fclose);
}

////////////////////////////////////////////////////////////////////

vector<string> H5Core::openFiles()
{
    Lock lock;
    std::set<string> names;
    forEachOpenObject([&names](const string& name, unsigned int /*intent*/) { names.insert(name); });
    return vector<string>(names.begin(), names.end());
}

////////////////////////////////////////////////////////////////////

#else

namespace
{
    // throws the fatal error for an attempt to use HDF5 in a build without HDF5 support
    [[noreturn]] void throwUnavailable()
    {
        throw FATALERROR("HDF5 files are not supported because this version of the code was built without HDF5\n"
                         "(enable the BUILD_HDF5 build option and rebuild the code)");
    }
}

////////////////////////////////////////////////////////////////////

bool H5Core::available()
{
    return false;
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::openGroup(const string& /*filePath*/, const string& /*groupPath*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

bool H5Core::hasGroup(const string& /*filePath*/, const string& /*groupPath*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::openGroup(std::int64_t /*group*/, const string& /*name*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

vector<string> H5Core::groupNames(std::int64_t /*group*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

vector<string> H5Core::datasetNames(std::int64_t /*group*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

bool H5Core::hasDataset(std::int64_t /*group*/, const string& /*name*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::openDataset(std::int64_t /*group*/, const string& /*name*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

void H5Core::close(std::int64_t /*object*/) noexcept {}

////////////////////////////////////////////////////////////////////

bool H5Core::hasAttribute(std::int64_t /*object*/, const string& /*name*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

bool H5Core::boolAttribute(std::int64_t /*object*/, const string& /*name*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

int H5Core::intAttribute(std::int64_t /*object*/, const string& /*name*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

double H5Core::doubleAttribute(std::int64_t /*object*/, const string& /*name*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

string H5Core::stringAttribute(std::int64_t /*object*/, const string& /*name*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

vector<size_t> H5Core::shape(std::int64_t /*dataset*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

vector<char> H5Core::readBool(std::int64_t /*dataset*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

vector<int> H5Core::readInt(std::int64_t /*dataset*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

Array H5Core::readDouble(std::int64_t /*dataset*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::createGroup(const string& /*filePath*/, const string& /*groupPath*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

void H5Core::setBoolAttribute(std::int64_t /*object*/, const string& /*name*/, bool /*value*/,
                              const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

void H5Core::setIntAttribute(std::int64_t /*object*/, const string& /*name*/, int /*value*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

void H5Core::setDoubleAttribute(std::int64_t /*object*/, const string& /*name*/, double /*value*/,
                                const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

void H5Core::setStringAttribute(std::int64_t /*object*/, const string& /*name*/, const string& /*value*/,
                                const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::createBoolDataset(std::int64_t /*group*/, const string& /*name*/, const vector<size_t>& /*shape*/,
                                       const char* /*values*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::createIntDataset(std::int64_t /*group*/, const string& /*name*/, const vector<size_t>& /*shape*/,
                                      const int* /*values*/, const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::createDoubleDataset(std::int64_t /*group*/, const string& /*name*/,
                                         const vector<size_t>& /*shape*/, const double* /*values*/,
                                         const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

std::int64_t H5Core::createTextDataset(std::int64_t /*group*/, const string& /*name*/, const string& /*text*/,
                                       const string& /*context*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

void H5Core::verifyWritable(const string& /*filePath*/)
{
    throwUnavailable();
}

////////////////////////////////////////////////////////////////////

vector<string> H5Core::openFiles()
{
    return vector<string>();
}

////////////////////////////////////////////////////////////////////

#endif
