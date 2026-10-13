/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StoredTableImpl.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "H5Lib.hpp"
#include "Log.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <cstring>

////////////////////////////////////////////////////////////////////

namespace
{
    // the alternate interpretations for 8-byte items in the stored table format
    union StabItem
    {
        double doubleType;
        size_t sizeType;
        char stringType[8];
    };
    const size_t itemSize = sizeof(StabItem);

    static_assert((sizeof(size_t) == 8) & (sizeof(double) == 8) & (itemSize == 8),
                  "Cannot properly declare union for items in stored table format");
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the name part of the given specification string, which has the form "name(unit)"
    string nameOf(const string& specification)
    {
        auto index = specification.find('(');
        if (index == string::npos || index == 0 || index > 8)
            throw FATALERROR("Invalid stored table name specification: " + specification);
        return specification.substr(0, index);
    }

    // returns the unit part of the given specification string, which has the form "name(unit)"
    string unitOf(const string& specification)
    {
        auto size = specification.size();
        auto index = specification.find('(');
        if (index == string::npos || size - index <= 2 || size - index > 10 || specification[size - 1] != ')')
            throw FATALERROR("Invalid stored table unit specification: " + specification);
        return specification.substr(index + 1, size - index - 2);
    }

    // returns true if the name part of the given specification string matches the given 8-byte item
    bool matchesName(const string& specification, const StabItem* nameItem)
    {
        string nameString = StringUtils::padRight(nameOf(specification), itemSize);
        return memcmp(nameString.c_str(), nameItem->stringType, itemSize) == 0;
    }

    // returns true if the unit part of the given specification string matches the given 8-byte item
    bool matchesUnit(const string& specification, const StabItem* unitItem)
    {
        string unitString = StringUtils::padRight(unitOf(specification), itemSize);
        return memcmp(unitString.c_str(), unitItem->stringType, itemSize) == 0;
    }
}

namespace
{
    // parses a stored table starting at the given memory location, verifying that it has the given
    // number of axes and matches the given axes and quantity specifications, and storing the results
    // in the given output parameters; the given label is used only for error messages
    void parse(const void* start, size_t numAxes, const string& axes, const string& quantity, const string& label,
               const double** axBeg, const double** qtyBeg, size_t* axLen, size_t* qtyStep, bool* axLog, bool* qtyLog)
    {
        const StabItem* currentItem = static_cast<const StabItem*>(start);

        // verify the name tag and the Endianness tag
        if (memcmp("SKIRT X\n", currentItem++->stringType, itemSize) || currentItem++->sizeType != 0x010203040A0BFEFF)
            throw FATALERROR("File does not have stored table format: " + label);

        // verify the number of axes
        if (currentItem++->sizeType != numAxes)
            throw FATALERROR("Number of axes in stored table file does not match: " + label);

        // split the axes specification into a list
        auto axesSpecs = StringUtils::split(axes, ",");
        if (axesSpecs.size() != numAxes)
            throw FATALERROR("Number of axes in stored table axes specification does not match: " + axes);

        // verify the axes names
        for (size_t i = 0; i < numAxes; ++i)
        {
            if (!matchesName(axesSpecs[i], currentItem++))
                throw FATALERROR("Axis " + std::to_string(i) + " does not have name " + axesSpecs[i]
                                 + " in stored table " + label);
        }

        // verify the axes units
        for (size_t i = 0; i < numAxes; ++i)
        {
            if (!matchesUnit(axesSpecs[i], currentItem++))
                throw FATALERROR("Axis " + std::to_string(i) + " does not have unit " + axesSpecs[i]
                                 + " in stored table " + label);
        }

        // store the interpolation scale for each axis
        for (size_t i = 0; i < numAxes; ++i)
        {
            *axLog++ = (memcmp("log     ", currentItem++->stringType, itemSize) == 0);
        }

        // store the grid length and a pointer to the first value for each axis
        for (size_t i = 0; i < numAxes; ++i)
        {
            size_t length = currentItem++->sizeType;
            *axLen++ = length;
            *axBeg++ = &currentItem->doubleType;
            currentItem += length;
        }

        // look for the appropriate quantity
        size_t numQties = currentItem++->sizeType;
        size_t qtyIndex = numQties;
        for (size_t i = 0; i < numQties; ++i)
        {
            if (matchesName(quantity, currentItem++)) qtyIndex = i;
        }
        if (qtyIndex == numQties)
            throw FATALERROR("Tabulated quantity " + quantity + " is not in stored table " + label);

        // verify the corresponding unit
        for (size_t i = 0; i < numQties; ++i)
        {
            if (i == qtyIndex)
            {
                if (!matchesUnit(quantity, currentItem++))
                    throw FATALERROR("Tabulated quantity does not have unit " + quantity + " in stored table " + label);
            }
            else
                currentItem++;
        }

        // store the corresponding interpolation scale
        for (size_t i = 0; i < numQties; ++i)
        {
            if (i == qtyIndex)
                *qtyLog = (memcmp("log     ", currentItem++->stringType, itemSize) == 0);
            else
                currentItem++;
        }

        // calculate and store the pointer to the first quantity value, and store the number of quantities
        *qtyBeg = &currentItem->doubleType + qtyIndex;
        *qtyStep = numQties;
    }

    // opens the dataset that represents the axis or quantity with the given specification in the given stored table
    // bundle, verifies that its unit matches the specification, and returns the dataset; a vertical bar in the name
    // of a dataset stands for a forward slash in the name of the axis or quantity; the last argument describes the
    // role of the dataset for use in error messages
    H5DatasetR openTableDataset(const H5BundleR& bundle, const string& specification, const string& label,
                                const string& role)
    {
        string name = StringUtils::replace(nameOf(specification), "/", "|");
        if (!bundle.hasDataset(name))
            throw FATALERROR(role + " " + specification + " is not in stored table bundle " + label);
        H5DatasetR dataset = bundle.getDataset(name);
        if (!dataset.hasAttribute("unit") || dataset.getStringAttribute("unit") != unitOf(specification))
            throw FATALERROR(role + " does not have unit " + specification + " in stored table bundle " + label);
        if (!dataset.hasAttribute("log"))
            throw FATALERROR(role + " " + specification + " has no 'log' attribute in stored table bundle " + label);
        return dataset;
    }

    // loads a stored table from the bundle with the given location, verifying that it has the given number of axes
    // and matches the given axes and quantity specifications; the grid points of the axes and the values of the
    // quantity are copied to the given storage, in the layout of a stored table file with a single quantity, and
    // the other output parameters are set to refer to that storage
    void load(const string& bundlePath, size_t numAxes, const string& axes, const string& quantity,
              vector<double>& storage, const double** axBeg, const double** qtyBeg, size_t* axLen, size_t* qtyStep,
              bool* axLog, bool* qtyLog)
    {
        H5BundleR bundle = H5Lib::openBundle(bundlePath);

        // split the axes specification into a list
        auto axesSpecs = StringUtils::split(axes, ",");
        if (axesSpecs.size() != numAxes)
            throw FATALERROR("Number of axes in stored table axes specification does not match: " + axes);

        // read the grid points for each axis, verifying the position of the axis and its unit
        vector<Array> grids(numAxes);
        size_t numPoints = 0;
        size_t numValues = 1;
        for (size_t i = 0; i < numAxes; ++i)
        {
            string role = "Axis " + std::to_string(i);
            H5DatasetR dataset = openTableDataset(bundle, axesSpecs[i], bundlePath, role);
            if (!dataset.hasAttribute("axis") || dataset.getIntAttribute("axis") != static_cast<int>(i))
                throw FATALERROR(role + " does not have name " + axesSpecs[i] + " in stored table bundle "
                                 + bundlePath);
            if (dataset.shape().size() != 1 || dataset.shape()[0] == 0)
                throw FATALERROR(role + " is not a nonempty one-dimensional dataset in stored table bundle "
                                 + bundlePath);
            axLog[i] = dataset.getBoolAttribute("log");
            grids[i] = dataset.readDouble();
            axLen[i] = grids[i].size();
            numPoints += grids[i].size();
            numValues *= grids[i].size();
        }

        // read the values of the quantity, verifying its unit and its shape; a dataset with an axis attribute
        // represents an axis rather than a quantity
        H5DatasetR dataset = openTableDataset(bundle, quantity, bundlePath, "Tabulated quantity");
        if (dataset.hasAttribute("axis"))
            throw FATALERROR("Tabulated quantity " + quantity + " is not in stored table bundle " + bundlePath);
        vector<size_t> shape = dataset.shape();
        if (shape.size() != numAxes)
            throw FATALERROR("Number of axes in stored table bundle does not match: " + bundlePath);
        for (size_t i = 0; i < numAxes; ++i)
        {
            if (shape[i] != axLen[i])
                throw FATALERROR("Tabulated quantity " + quantity + " does not match the length of axis "
                                 + std::to_string(i) + " in stored table bundle " + bundlePath);
        }
        *qtyLog = dataset.getBoolAttribute("log");
        Array values = dataset.readDouble();

        // copy the grid points for each axis to the storage, one axis after the other
        storage.resize(numPoints + numValues);
        double* target = storage.data();
        for (size_t i = 0; i < numAxes; ++i)
        {
            axBeg[i] = target;
            target = std::copy(begin(grids[i]), end(grids[i]), target);
        }

        // copy the values to the storage, reordering them from the layout of the dataset, where the index for the
        // last axis varies most rapidly, to the layout of a stored table, where the index for the first axis does
        *qtyBeg = target;
        *qtyStep = 1;
        vector<size_t> indices(numAxes);  // the current index for each axis, in the layout of the dataset
        for (size_t source = 0; source != numValues; ++source)
        {
            // calculate the position of this value in the layout of a stored table
            size_t position = indices[numAxes - 1];
            for (size_t k = numAxes - 2; k < numAxes; --k) position = position * axLen[k] + indices[k];
            target[position] = values[source];

            // advance the indices to the next value in the layout of the dataset
            for (size_t k = numAxes - 1; k < numAxes; --k)
            {
                if (++indices[k] != axLen[k]) break;
                indices[k] = 0;
            }
        }
    }

    // logs that the given item opened the stored table identified by the given label, unless the same
    // thread already successfully opened the same label for the same item
    void logOpened(const SimulationItem* item, const string& label)
    {
        thread_local const SimulationItem* previousItem = nullptr;
        thread_local string previousLabel;
        if (item != previousItem || label != previousLabel)
        {
            item->find<Log>()->info(item->type() + " opened stored table " + label);
            previousItem = item;
            previousLabel = label;
        }
    }
}

////////////////////////////////////////////////////////////////////

void StoredTable_Impl::open(size_t numAxes, const SimulationItem* item, string filename, bool resource,
                            const string& axes, const string& quantity, string& filePath, vector<double>& storage,
                            const double** axBeg, const double** qtyBeg, size_t* axLen, size_t* qtyStep, bool* axLog,
                            bool* qtyLog)
{
    // add the mandatory filename extension if needed
    if (!StringUtils::endsWith(filename, ".stab")) filename += ".stab";

    string mapPath;
    if (resource)
    {
        // retrieve the full path for the resource; the function throws a fatal error if the resource cannot be found
        mapPath = FilePaths::resource(filename);
    }
    else
    {
        // determine the candidate locations for the input file: a regular file and a bundle in an HDF5 file
        auto [candidatePath, bundlePath] = item->find<FilePaths>()->input(filename);

        // load the table from the bundle if the regular file does not exist and there is a bundle candidate
        if (!bundlePath.empty() && (candidatePath.empty() || !System::isFile(candidatePath)))
        {
            load(bundlePath, numAxes, axes, quantity, storage, axBeg, qtyBeg, axLen, qtyStep, axLog, qtyLog);
            logOpened(item, bundlePath);
            return;
        }
        mapPath = candidatePath;
    }

    // acquire a memory map for the file; the function returns zeros if the memory map cannot be created
    auto [start, size] = System::acquireMemoryMap(mapPath);
    if (!start) throw FATALERROR("Cannot acquire memory map for file: " + mapPath);
    filePath = mapPath;

    parse(start, numAxes, axes, quantity, filePath, axBeg, qtyBeg, axLen, qtyStep, axLog, qtyLog);
    logOpened(item, filePath);
}

////////////////////////////////////////////////////////////////////

void StoredTable_Impl::openAt(size_t numAxes, const string& filePath, size_t byteOffset, const string& label,
                              const string& axes, const string& quantity, const double** axBeg, const double** qtyBeg,
                              size_t* axLen, size_t* qtyStep, bool* axLog, bool* qtyLog)
{
    // acquire a memory map for the file; the function returns zeros if the memory map cannot be created
    auto [start, size] = System::acquireMemoryMap(filePath);
    if (!start) throw FATALERROR("Cannot acquire memory map for file: " + filePath);

    // unlike open(), this function does not log anything on success: the dictionary that spawned
    // this stored table already logged a single message when it was opened
    parse(static_cast<const char*>(start) + byteOffset, numAxes, axes, quantity, label, axBeg, qtyBeg, axLen, qtyStep,
          axLog, qtyLog);
}

////////////////////////////////////////////////////////////////////

void StoredTable_Impl::close(const std::string& filePath)
{
    if (!filePath.empty()) System::releaseMemoryMap(filePath);
}

////////////////////////////////////////////////////////////////////
