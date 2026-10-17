/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ColumnInFile.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "H5Lib.hpp"
#include "InputColumns.hpp"
#include "Log.hpp"
#include "SimulationItem.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"
#include <exception>

////////////////////////////////////////////////////////////////////

ColumnInFile::ColumnInFile(const SimulationItem* item, const string& filename, const string& description)
    : _item(item), _filename(filename), _description(description)
{}

////////////////////////////////////////////////////////////////////

void ColumnInFile::close()
{
    if (_isClosed) return;
    _isClosed = true;

    // close the regular file, which logs its own "done" message
    if (_text)
    {
        _text->close();
    }
    // or release the data read from the bundle, and log a "done" message except after an exception has been thrown
    else if (_isOpen)
    {
        _columns.clear();
        if (std::uncaught_exceptions() == 0) _log->info("Done reading");
    }
}

////////////////////////////////////////////////////////////////////

ColumnInFile::~ColumnInFile()
{
    close();
}

////////////////////////////////////////////////////////////////////

void ColumnInFile::useColumns(const string& columns)
{
    if (_isOpen || _isClosed) throw FATALERROR("Column remapping was requested after reading from the input file");
    _declarations.push_back({true, columns, string(), string()});
}

////////////////////////////////////////////////////////////////////

void ColumnInFile::addColumn(const string& description, const string& quantity, const string& defaultUnit)
{
    if (_isOpen || _isClosed) throw FATALERROR("Program columns were declared after reading from the input file");
    _declarations.push_back({false, description, quantity, defaultUnit});
}

////////////////////////////////////////////////////////////////////

void ColumnInFile::open()
{
    if (_isOpen || _isClosed) return;
    _isOpen = true;

    // the binary stored columns format, recognized by its filename extension, is no longer supported
    if (StringUtils::endsWith(_filename, ".scol"))
        throw FATALERROR("The stored columns format of input file " + _filename
                         + " is no longer supported; convert the file to a column text file or to a bundle in an "
                           "HDF5 file with PTS, and refer to the result in the ski file");

    // determine the candidate locations for the input file: a regular file and a bundle in an HDF5 file
    auto [filePath, bundlePath] = _item->find<FilePaths>()->input(_filename);

    // read from the regular file if it exists, or if there is no bundle candidate
    if (bundlePath.empty() || (!filePath.empty() && System::isFile(filePath)))
    {
        _text = std::make_unique<TextInFile>(_item, filePath, _description);
        for (const Declaration& declaration : _declarations)
        {
            if (declaration.isUseColumns)
                _text->useColumns(declaration.first);
            else
                _text->addColumn(declaration.first, declaration.quantity, declaration.unit);
        }
    }
    // otherwise read from the bundle; if there is no bundle either, report both locations
    else
    {
        _item->find<FilePaths>()->verifyBundleInput(bundlePath);
        if (!filePath.empty() && !H5Lib::hasBundle(bundlePath))
            throw FATALERROR("Could not find the " + _description + " input file " + _filename
                             + ": there is no regular file " + filePath + " and no bundle " + bundlePath);
        _bundlePath = bundlePath;
        openBundle();
    }
}

////////////////////////////////////////////////////////////////////

void ColumnInFile::openBundle()
{
    // open the bundle
    _log = _item->find<Log>();
    H5BundleR bundle = H5Lib::openBundle(_bundlePath);
    _log->info(_item->typeAndName() + " reads " + _description + " from bundle " + _bundlePath + "...");

    // allow a subclass to read extra datasets
    openedBundle(bundle);

    // locate the dataset for each physical column: a dataset with a "column" attribute holding the index of the column
    vector<string> datasetNames;
    for (const string& name : bundle.getDatasetNames())
    {
        H5DatasetR dataset = bundle.getDataset(name);
        if (dataset.hasAttribute("column"))
        {
            int column = dataset.getIntAttribute("column");
            if (column < 0)
                throw FATALERROR("Dataset '" + name + "' in bundle " + _bundlePath + " has a negative column index");
            size_t index = static_cast<size_t>(column);
            if (datasetNames.size() <= index) datasetNames.resize(index + 1);
            if (!datasetNames[index].empty())
                throw FATALERROR("Datasets '" + datasetNames[index] + "' and '" + name + "' in bundle " + _bundlePath
                                 + " have the same column index " + std::to_string(column));
            datasetNames[index] = name;
        }
    }
    for (size_t index = 0; index != datasetNames.size(); ++index)
    {
        if (datasetNames[index].empty())
            throw FATALERROR("There is no dataset with column index " + std::to_string(index) + " in bundle "
                             + _bundlePath + "; the column indices must be consecutive and start at zero");
    }

    // register the physical columns; a vertical bar in the name of a dataset stands for a forward slash
    InputColumns columns(_item->find<Units>(), _log, "bundle");
    for (const string& name : datasetNames)
    {
        H5DatasetR dataset = bundle.getDataset(name);
        string title = StringUtils::squeeze(StringUtils::replace(name, "|", "/"));
        if (dataset.hasAttribute("unit"))
            columns.addFileColumn(title, StringUtils::squeeze(dataset.getStringAttribute("unit")));
        else
            columns.addFileColumn(title);
    }

    // register the logical columns expected by the program, which determines the columns to be read
    for (const Declaration& declaration : _declarations)
    {
        if (declaration.isUseColumns)
            columns.useColumns(declaration.first);
        else
        {
            if (datasetNames.empty())
                throw FATALERROR("Bundle " + _bundlePath + " does not hold any datasets with a column index");
            columns.addColumn(declaration.first, declaration.quantity, declaration.unit);
        }
    }

    // read the dataset for each physical column that is needed, in the order of these columns, and convert its values
    // to internal units; the order matters because the conversion may depend on a preceding wavelength column
    _numLogCols = columns.numLogicalColumns();
    _columns.resize(_numLogCols);
    const vector<size_t>& logicalIndices = columns.logicalIndices();
    for (size_t index = 0; index != logicalIndices.size(); ++index)
    {
        size_t logicalIndex = logicalIndices[index];
        if (logicalIndex != InputColumns::noIndex)
        {
            const string& name = datasetNames[index];
            H5DatasetR dataset = bundle.getDataset(name);
            if (dataset.shape().size() != 1)
                throw FATALERROR("Dataset '" + name + "' in bundle " + _bundlePath + " is not one-dimensional");
            Array values = dataset.readDouble();

            // all columns must have the same number of rows
            if (!_hasNumRows)
            {
                _numRows = values.size();
                _hasNumRows = true;
            }
            else if (values.size() != _numRows)
            {
                throw FATALERROR("Dataset '" + name + "' in bundle " + _bundlePath + " has "
                                 + std::to_string(values.size()) + " rows rather than " + std::to_string(_numRows));
            }

            columns.convertColumn(logicalIndex, values, _columns);
            _columns[logicalIndex] = std::move(values);
        }
    }

    // if no columns were read, obtain the number of rows from the first column in the bundle, if there is one
    if (!_hasNumRows && !datasetNames.empty())
    {
        vector<size_t> shape = bundle.getDataset(datasetNames[0]).shape();
        if (shape.size() != 1)
            throw FATALERROR("Dataset '" + datasetNames[0] + "' in bundle " + _bundlePath + " is not one-dimensional");
        _numRows = shape[0];
    }

    // provide zero values for the virtual columns, which do not correspond to a physical column
    for (size_t logicalIndex = 0; logicalIndex != _numLogCols; ++logicalIndex)
    {
        if (std::find(logicalIndices.begin(), logicalIndices.end(), logicalIndex) == logicalIndices.end())
            _columns[logicalIndex].resize(_numRows);
    }
}

////////////////////////////////////////////////////////////////////

void ColumnInFile::openedBundle(const H5BundleR& /*bundle*/) {}

////////////////////////////////////////////////////////////////////

void ColumnInFile::setNumRows(size_t numRows)
{
    _numRows = numRows;
    _hasNumRows = true;
}

////////////////////////////////////////////////////////////////////

bool ColumnInFile::readRow(Array& values)
{
    open();
    if (_text) return _text->readRow(values);

    // read the next row from the columns loaded from the bundle, if there is one
    if (_isClosed || _nextRow == _numRows) return false;
    if (values.size() != _numLogCols) values.resize(_numLogCols);
    for (size_t index = 0; index != _numLogCols; ++index) values[index] = _columns[index][_nextRow];
    _nextRow++;
    return true;
}

////////////////////////////////////////////////////////////////////

vector<Array> ColumnInFile::readAllRows()
{
    open();
    if (_text) return _text->readAllRows();

    vector<Array> rows;
    if (!_isClosed) rows.reserve(_numRows - _nextRow);
    while (true)
    {
        rows.emplace_back();        // add a default-constructed array to the vector
        if (!readRow(rows.back()))  // read next row's values into that array
        {
            rows.pop_back();  // at the end, remove the extraneous array
            break;
        }
    }
    return rows;
}

////////////////////////////////////////////////////////////////////

vector<Array> ColumnInFile::readAllColumns()
{
    open();
    if (_text) return _text->readAllColumns();

    // hand over the columns loaded from the bundle, or the remaining portion if some rows were already read
    vector<Array> columns(_numLogCols);
    if (!_isClosed)
    {
        if (_nextRow == 0)
        {
            columns = std::move(_columns);
            _columns.clear();
        }
        else
        {
            for (size_t index = 0; index != _numLogCols; ++index)
                columns[index] = Array(_columns[index][std::slice(_nextRow, _numRows - _nextRow, 1)]);
        }
        _nextRow = _numRows;
    }
    return columns;
}

////////////////////////////////////////////////////////////////////
