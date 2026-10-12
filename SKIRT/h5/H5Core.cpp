/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "H5Core.hpp"
#include "FatalError.hpp"

#ifdef BUILD_HDF5
#    include "System.hpp"
#    include <hdf5.h>
#    include <mutex>
#endif

////////////////////////////////////////////////////////////////////

#ifdef BUILD_HDF5

namespace
{
    // the handles in the interface of this namespace must be able to hold the identifiers used by the HDF5 library
    static_assert(sizeof(hid_t) <= sizeof(std::int64_t), "HDF5 identifiers do not fit in a 64-bit integer");

    // the mutex that protects all calls into the HDF5 library
    std::mutex _mutex;

    // An instance of this class locks the mutex for its lifetime. It also stops the HDF5 library from printing its
    // own error messages for the calling thread, because failures are reported through fatal errors instead.
    class Lock
    {
        std::unique_lock<std::mutex> _lock;

    public:
        Lock() : _lock(_mutex) { H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr); }
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

    // open the file for reading; on a file system that does not support file locking, proceed without locking
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

    // open the group; the library keeps the file open as long as the group is open
    Id group(H5Gopen2(file.get(), ("/" + groupPath).c_str(), H5P_DEFAULT), H5Gclose);
    if (!group.valid()) throw FATALERROR("HDF5 file " + filePath + " has no group '" + groupPath + "'");
    return group.release();
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

vector<string> H5Core::groupNames(std::int64_t group, const string& context)
{
    Lock lock;

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

        // include the name only if it refers to a group
        Id object(H5Oopen(group, name.c_str(), H5P_DEFAULT), H5Oclose);
        if (object.valid() && H5Iget_type(object.get()) == H5I_GROUP) names.push_back(name);
    }
    return names;
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

#endif
