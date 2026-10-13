/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InputColumns.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // Error return values for the functions below, in addition to InputColumns::noIndex
    const int ERROR_NO_EXPON = 999;
    const size_t ERROR_AM_INDEX = 999998;

    // This function returns the wavelength exponent needed to convert a per wavelength/frequency/energy
    // quantity to internal (per wavelength) style, given the input units, or the error value if the
    // given units are not supported by any of the relevant quantities.
    int waveExponentForSpecificQuantity(Units* unitSystem, const string& unitString)
    {
        // a list of known per wavelength / per frequency quantities and the corresponding exponents
        static const vector<string> specificQuantities(
            {"wavelengthmonluminosity", "wavelengthfluxdensity", "wavelengthsurfacebrightness", "neutralmonluminosity",
             "neutralfluxdensity", "neutralsurfacebrightness", "frequencymonluminosity", "frequencyfluxdensity",
             "frequencysurfacebrightness", "energymonluminosity", "energyfluxdensity", "energysurfacebrightness"});
        static const vector<int> specificExponents({0, 0, 0, -1, -1, -1, -2, -2, -2, -3, -3, -3});

        // loop over the list
        for (size_t q = 0; q != specificQuantities.size(); ++q)
        {
            // if this quantity supports the given unit, return the corresponding exponent
            if (unitSystem->has(specificQuantities[q], unitString)) return specificExponents[q];
        }
        return ERROR_NO_EXPON;
    }
}

////////////////////////////////////////////////////////////////////

InputColumns::InputColumns(Units* units, Log* log, string source) : _units(units), _log(log), _source(std::move(source))
{}

////////////////////////////////////////////////////////////////////

void InputColumns::addFileColumn(const string& name, const string& unit)
{
    // add a default-constructed ColumnInfo record to the list
    _colv.emplace_back();

    // remember the description and the units specified in the file
    _colv.back().physColIndex = _colv.size();  // one-based column index
    _colv.back().unit = unit;
    _colv.back().title = name;
    _hasFileInfo = true;
}

////////////////////////////////////////////////////////////////////

void InputColumns::addFileColumn(const string& name)
{
    addFileColumn(name, string());
    _colv.back().hasFileUnit = false;
}

////////////////////////////////////////////////////////////////////

size_t InputColumns::indexForName(const std::string& name) const
{
    size_t result = noIndex;
    size_t index = 0;
    for (const ColumnInfo& col : _colv)
    {
        if (col.title == name)
        {
            if (result != noIndex) return ERROR_AM_INDEX;
            result = index;
        }
        index++;
    }
    return result;
}

////////////////////////////////////////////////////////////////////

size_t InputColumns::waveIndexForSpecificQuantity() const
{
    for (size_t index = 0; index < _numLogCols; ++index)
    {
        if (_colv[index].description == "wavelength")
        {
            if (_colv[index].physColIndex < _colv[_numLogCols - 1].physColIndex)
                return index;
            else
                break;
        }
    }
    return noIndex;
}

////////////////////////////////////////////////////////////////////

void InputColumns::useColumns(string columns)
{
    // empty columns string behaves as if we were never called at all
    columns = StringUtils::squeeze(columns);
    if (columns.empty()) return;

    // verify that program columns have not yet been added
    if (_hasProgInfo) throw FATALERROR("Program columns were declared before requesting column remapping");

    // verify that file contains column info
    if (!_hasFileInfo) throw FATALERROR("Requesting logical columns but there is no column info in " + _source);

    // automatic column assignment happens in addColumn() because we don't have the required information here
    if (columns == "*" || columns == "*0")
    {
        _doAutoCols = true;
        if (columns == "*0") _allowZeroCols = true;
        return;
    }

    // establish the logical column info list
    vector<ColumnInfo> newcolv;
    for (const string& name : StringUtils::split(columns, ","))
    {
        string sname = StringUtils::squeeze(name);

        // handle virtual zero column
        if (sname == "0")
        {
            // add a default-constructed ColumnInfo record; physColIndex remains zero to indicate virtual zero column
            newcolv.emplace_back();
            _haveZeroCols = true;
        }
        // handle regular physical column
        else
        {
            size_t index = indexForName(sname);
            if (index == noIndex)
                throw FATALERROR("No column description in " + _source + " matches logical name '" + sname + "'");
            if (index == ERROR_AM_INDEX)
                throw FATALERROR("Multiple column descriptions in " + _source + " match logical name '" + sname + "'");

            // copy the appropriate ColumnInfo record
            newcolv.emplace_back(_colv[index]);
        }
    }

    // replace the column info list
    _colv = newcolv;
}

////////////////////////////////////////////////////////////////////

void InputColumns::addColumn(const string& description, const string& quantity, const string& defaultUnit)
{
    _hasProgInfo = true;

    // in case of automatic column assignment, find the intended column and move it to the correct place
    if (_doAutoCols)
    {
        string sdescription = StringUtils::squeeze(description);
        size_t index = indexForName(sdescription);
        if (index == noIndex)
        {
            if (_allowZeroCols)
            {
                // insert a default-constructed record; physColIndex remains zero to indicate virtual zero column
                _colv.emplace(_colv.begin() + _numLogCols);
                _haveZeroCols = true;
            }
            else
                throw FATALERROR("No column description in " + _source + " matches '" + sdescription + "'");
        }
        else if (index == ERROR_AM_INDEX)
            throw FATALERROR("Multiple column descriptions in " + _source + " match '" + sdescription + "'");
        else if (index < _numLogCols)
            throw FATALERROR("Program requests multiple columns with same description '" + sdescription + "'");
        else if (index > _numLogCols)
            std::swap(_colv[index], _colv[_numLogCols]);
    }
    // if the file has no column info at all, add a default record for this column
    else if (!_hasFileInfo)
    {
        _colv.emplace_back();
        _colv.back().physColIndex = _numLogCols + 1;
        _colv.back().unit = defaultUnit;
    }
    // otherwise verify that there is a column specification to match this program column
    else
    {
        if (_numLogCols + 1 > _colv.size())
            throw FATALERROR("No column info in " + _source + " for column " + std::to_string(_numLogCols + 1));
    }

    // get a writable reference to the column record being handled, and increment the program column index
    ColumnInfo& col = _colv[_numLogCols++];

    // store the programmatically provided information in the record;
    // the unit is already stored, unless the file has information on the column but does not specify its unit
    col.description = description;
    col.quantity = quantity;
    if (!col.hasFileUnit) col.unit = defaultUnit;

    // if this is not a virtual zero column, verify units and add physical to logical column mapping
    if (col.physColIndex)
    {
        // verify units and determine conversion factor for this column
        if (col.quantity.empty())  // dimensionless quantity
        {
            if (!col.unit.empty() && col.unit != "1")
                throw FATALERROR("Invalid units (" + col.unit + ") for dimensionless quantity in column "
                                 + std::to_string(_numLogCols));
            col.unit = "1";
        }
        else if (col.quantity == "specific")  // arbitrarily scaled value per wavelength or per frequency
        {
            col.waveExponent = waveExponentForSpecificQuantity(_units, col.unit);
            if (col.waveExponent == ERROR_NO_EXPON)
                throw FATALERROR("Invalid units (" + col.unit + ") for specific quantity '" + col.quantity
                                 + "' in column " + std::to_string(_numLogCols));
            if (col.waveExponent)
            {
                col.waveIndex = waveIndexForSpecificQuantity();
                if (col.waveIndex == noIndex)
                    throw FATALERROR("No preceding wavelength column for specific quantity '" + col.quantity
                                     + "' in column " + std::to_string(_numLogCols));
            }
        }
        else
        {
            if (!_units->has(col.quantity, col.unit))
                throw FATALERROR("Invalid units (" + col.unit + ") for quantity '" + col.quantity + "' in column "
                                 + std::to_string(_numLogCols));
            double offset;  // all SKIRT units have a zero offset
            std::tie(col.convFactor, col.convPower, offset) = _units->def(col.quantity, col.unit);
        }

        // add the physical to logical column mapping for this column
        if (_logColIndices.size() < col.physColIndex) _logColIndices.resize(col.physColIndex, noIndex);
        if (_logColIndices[col.physColIndex - 1] != noIndex)
            throw FATALERROR("Multiple logical columns (" + std::to_string(_logColIndices[col.physColIndex - 1] + 1)
                             + "," + std::to_string(_numLogCols) + ") map to the same physical column ("
                             + std::to_string(col.physColIndex) + ")");
        _logColIndices[col.physColIndex - 1] = _numLogCols - 1;
    }

    // log column information if requested
    if (_log)
    {
        string message = "  Column " + std::to_string(_numLogCols) + ": " + col.description;

        // for a regular physical column
        if (col.physColIndex)
        {
            if (!col.title.empty())
            {
                message += " <-- ";
                if (col.physColIndex != _numLogCols) message += "column " + std::to_string(col.physColIndex) + ": ";
                message += col.title;
            }
            message += " (" + col.unit + ")";
        }
        // for a virtual zero column
        else
        {
            message += " <-- 0";
        }
        _log->info(message);
    }
}

////////////////////////////////////////////////////////////////////

void InputColumns::convertColumn(size_t logicalIndex, Array& values, const vector<Array>& columns) const
{
    const ColumnInfo& col = _colv[logicalIndex];
    const Array* wavelengths = col.waveExponent ? &columns[col.waveIndex] : nullptr;

    size_t numRows = values.size();
    for (size_t r = 0; r != numRows; ++r)
    {
        double value = values[r];
        if (col.convPower != 1.) value = pow(value, col.convPower);
        values[r] = value * (wavelengths ? pow((*wavelengths)[r], col.waveExponent) : col.convFactor);
    }
}

////////////////////////////////////////////////////////////////////
