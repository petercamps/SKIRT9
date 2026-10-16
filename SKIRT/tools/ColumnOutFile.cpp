/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ColumnOutFile.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "H5Lib.hpp"
#include "Log.hpp"
#include "ProcessManager.hpp"
#include "SimulationItem.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "TextOutFile.hpp"
#include <exception>
#include <regex>

////////////////////////////////////////////////////////////////////

ColumnOutFile::ColumnOutFile(const SimulationItem* item, const string& filename, const string& description)
    : _item(item), _filename(filename), _description(description)
{}

////////////////////////////////////////////////////////////////////

void ColumnOutFile::close()
{
    if (_isClosed) return;

    // a file without rows still receives its header, or a bundle with empty datasets
    open();
    _isClosed = true;

    // close the regular file, which logs its own message
    if (_text)
    {
        _text->close();
        _text.reset();
    }

    // or write the bundle with the collected rows
    else if (!_bundlePath.empty())
    {
        string bundleDescription = _longDescription;
        vector<string> descriptions, units;
        for (const HeaderLine& line : _header)
        {
            if (line.isColumn)
            {
                descriptions.push_back(line.text);
                units.push_back(line.unit);
            }
            else
            {
                bundleDescription += (bundleDescription.empty() ? "" : "\n") + line.text;
            }
        }
        writeBundle(_item, _bundlePath, _description, bundleDescription, descriptions, units, _columns);
        _columns.clear();
    }
}

////////////////////////////////////////////////////////////////////

ColumnOutFile::~ColumnOutFile()
{
    // complete the output file if the client did not call close(), unless an exception is being propagated
    if (!_isClosed && std::uncaught_exceptions() == 0)
    {
        // a destructor must not throw, so log the message of a fatal error and ignore other exceptions
        try
        {
            try
            {
                close();
            }
            catch (FatalError& error)
            {
                auto log = _item->find<Log>();
                for (const string& line : error.message()) log->error(line);
            }
        }
        catch (...)
        {}
    }
}

////////////////////////////////////////////////////////////////////

void ColumnOutFile::setLongDescription(const string& longDescription)
{
    if (_isOpen) throw FATALERROR("The description of an output file was set after writing to the file");
    _longDescription = longDescription;
}

////////////////////////////////////////////////////////////////////

void ColumnOutFile::addColumn(const string& quantityDescription, string unitDescription, char format, int precision)
{
    if (_isOpen) throw FATALERROR("A column was added to an output file after writing to the file");
    _header.push_back({true, quantityDescription, unitDescription.empty() ? "1" : unitDescription, format, precision});
    _numColumns++;
}

////////////////////////////////////////////////////////////////////

void ColumnOutFile::addComment(const string& comment)
{
    if (_isOpen) throw FATALERROR("A comment was added to an output file after writing to the file");
    _header.push_back({false, comment, string(), 'e', 9});
}

////////////////////////////////////////////////////////////////////

void ColumnOutFile::writeRow(const vector<double>& values)
{
    writeRowPrivate(values.size(), values.data());
}

////////////////////////////////////////////////////////////////////

void ColumnOutFile::open()
{
    if (_isOpen) return;
    _isOpen = true;

    // only the root process writes the output file
    if (!ProcessManager::isRoot()) return;

    // determine the candidate locations for the output file: a regular file and a bundle in an HDF5 file
    auto bundlePath = _item->find<FilePaths>()->output(_filename + ".dat").second;

    // write a regular file through a TextOutFile instance, exactly as a client of that class would have written it
    if (bundlePath.empty())
    {
        _text = std::make_unique<TextOutFile>(_item, _filename, _description);
        if (!_longDescription.empty()) _text->writeLine("# " + _longDescription);
        for (const HeaderLine& line : _header)
        {
            if (line.isColumn)
                _text->addColumn(line.text, line.unit, line.format, line.precision);
            else
                _text->writeLine("# " + line.text);
        }
    }

    // or collect the rows for writing a bundle when the file is closed
    else
    {
        _bundlePath = bundlePath;
        _columns.resize(_numColumns);
    }
}

////////////////////////////////////////////////////////////////////

void ColumnOutFile::writeRowPrivate(size_t n, const double* values)
{
    if (n != _numColumns) throw FATALERROR("Number of values in row does not match the number of columns");
    if (_isClosed) return;
    open();

    if (_text)
        _text->writeRow(vector<double>(values, values + n));
    else if (!_bundlePath.empty())
        for (size_t i = 0; i != n; ++i) _columns[i].push_back(values[i]);
}

////////////////////////////////////////////////////////////////////

void ColumnOutFile::writeBundle(const SimulationItem* item, const string& bundlePath, const string& description,
                                const string& bundleDescription, const vector<string>& descriptions,
                                const vector<string>& units, vector<vector<double>>& columns)
{
    H5BundleW bundle = H5Lib::createBundle(bundlePath);
    if (!bundleDescription.empty()) bundle.setStringAttribute("description", bundleDescription);

    // derive a valid and unique dataset name from the description of each column
    vector<string> names;
    for (size_t index = 0; index != descriptions.size(); ++index)
    {
        string name = StringUtils::replace(StringUtils::squeeze(descriptions[index]), "/", "|");
        if (name.empty() || name == "." || name == "..") name = "column" + std::to_string(index);
        string unique = name;
        for (int count = 2; std::find(names.begin(), names.end(), unique) != names.end(); ++count)
            unique = name + " (" + std::to_string(count) + ")";
        names.push_back(unique);
    }

    // write a dataset for each column, releasing the memory for its values afterwards
    for (size_t index = 0; index != descriptions.size(); ++index)
    {
        Array values(columns[index].data(), columns[index].size());
        vector<double>().swap(columns[index]);
        H5DatasetW dataset = bundle.createDoubleDataset(names[index], values);
        dataset.setIntAttribute("column", static_cast<int>(index));
        dataset.setStringAttribute("unit", units[index]);
    }

    item->find<Log>()->info(item->typeAndName() + " wrote " + description + " to bundle " + bundlePath);
}

////////////////////////////////////////////////////////////////////

void ColumnOutFile::copyToBundle(const SimulationItem* item, const string& filename, const string& description)
{
    // only the root process writes output, and only if the output path includes an HDF5 file
    if (!ProcessManager::isRoot()) return;
    auto [filePath, bundlePath] = item->find<FilePaths>()->output(filename + ".dat");
    if (bundlePath.empty()) return;

    std::ifstream in = System::ifstream(filePath);
    if (!in) throw FATALERROR("Could not open the " + description + " text file " + filePath);

    // the syntax of a header line with column information, as for the TextInFile class
    static const std::regex syntax("#\\s*column\\s*(\\d*)\\s*:\\s*([^()]*)\\(\\s*([a-zA-Z0-9/]*)\\s*\\)\\s*",
                                   std::regex::icase);

    // read the header lines and the data lines
    string bundleDescription;
    vector<string> descriptions, units;
    vector<vector<double>> columns;
    string line;
    while (std::getline(in, line))
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        string trimmed = StringUtils::squeeze(line);
        if (trimmed.empty()) continue;

        // header line: column information, or a comment that becomes part of the description
        if (trimmed[0] == '#')
        {
            std::smatch matches;
            if (columns.empty() && std::regex_match(trimmed, matches, syntax))
            {
                descriptions.push_back(StringUtils::squeeze(matches[2].str()));
                units.push_back(matches[3].str());
            }
            else if (columns.empty())
            {
                string comment = StringUtils::squeeze(trimmed.substr(1));
                if (!comment.empty()) bundleDescription += (bundleDescription.empty() ? "" : "\n") + comment;
            }
            continue;
        }

        // data line: a value for each column
        if (columns.empty()) columns.resize(descriptions.size());
        auto tokens = StringUtils::split(trimmed, " ");
        if (tokens.size() < descriptions.size())
            throw FATALERROR("A data line in the " + description + " text file " + filePath + " has "
                             + std::to_string(tokens.size()) + " values rather than "
                             + std::to_string(descriptions.size()));
        for (size_t index = 0; index != descriptions.size(); ++index)
        {
            const char* begptr = tokens[index].c_str();
            char* endptr = nullptr;
            double value = std::strtod(begptr, &endptr);
            if (endptr == begptr || *endptr != '\0')
                throw FATALERROR("Value '" + tokens[index] + "' in the " + description + " text file " + filePath
                                 + " is not formatted as a floating point number");
            columns[index].push_back(value);
        }
    }
    if (descriptions.empty())
        throw FATALERROR("The " + description + " text file " + filePath + " has no column information");
    columns.resize(descriptions.size());

    writeBundle(item, bundlePath, description, bundleDescription, descriptions, units, columns);
}

////////////////////////////////////////////////////////////////////
