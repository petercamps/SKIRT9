/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TextOutFile.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "H5Lib.hpp"
#include "Log.hpp"
#include "ProcessManager.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "Units.hpp"
#include <exception>
#include <sstream>

////////////////////////////////////////////////////////////////////

TextOutFile::TextOutFile(const SimulationItem* item, const string& filename, const string& description,
                         const string& textType)
{
    // Only open the output file if this is the root process
    if (ProcessManager::isRoot())
    {
        // remember some pointers
        _log = item->find<Log>();
        _units = item->find<Units>();

        // determine the candidate locations for the output file: a regular file and a bundle in an HDF5 file
        auto [filepath, bundlePath] = item->find<FilePaths>()->output(filename + ".dat");

        // collect the text for a bundle if requested and possible
        if (!textType.empty() && !bundlePath.empty())
        {
            _bundlePath = bundlePath;
            _textType = textType;
            _message = item->typeAndName() + " wrote " + description + " to bundle " + bundlePath;
        }

        // otherwise open the regular file
        else
        {
            _out = System::ofstream(filepath);
            if (!_out) throw FATALERROR("Could not open the " + description + " output file " + filepath);
            _message = item->typeAndName() + " wrote " + description + " to " + filepath;
        }
    }
}

////////////////////////////////////////////////////////////////////

void TextOutFile::close()
{
    if (_out.is_open())
    {
        _out.close();

        // log success message, except if an exception has been thrown
        if (std::uncaught_exceptions() == 0) _log->info(_message);
    }
    else if (!_bundlePath.empty())
    {
        // write the bundle once, and release the collected text
        string bundlePath;
        string text;
        std::swap(bundlePath, _bundlePath);
        std::swap(text, _text);
        H5BundleW bundle = H5Lib::createBundle(bundlePath);
        bundle.setStringAttribute("type", _textType);
        bundle.createTextDataset("text", text);
        _log->info(_message);
    }
}

////////////////////////////////////////////////////////////////////

TextOutFile::~TextOutFile()
{
    // close a regular file in any case, but do not write a bundle while an exception is being propagated
    if (_out.is_open() || (!_bundlePath.empty() && std::uncaught_exceptions() == 0))
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
                for (const string& line : error.message()) _log->error(line);
            }
        }
        catch (...)
        {}
    }
}

////////////////////////////////////////////////////////////////////

string TextOutFile::copyToBundle(const SimulationItem* item, const string& name, const string& description,
                                 const string& textType)
{
    // only the root process writes output, and only if the output path includes an HDF5 file
    if (!ProcessManager::isRoot()) return string();
    auto [filePath, bundlePath] = item->find<FilePaths>()->output(name);
    if (bundlePath.empty()) return string();

    // read the complete file
    std::ifstream in = System::ifstream(filePath);
    if (!in) throw FATALERROR("Could not open the " + description + " " + filePath);
    std::ostringstream text;
    text << in.rdbuf();
    if (in.bad()) throw FATALERROR("Could not read the " + description + " " + filePath);

    // write the bundle
    H5BundleW bundle = H5Lib::createBundle(bundlePath);
    bundle.setStringAttribute("type", textType);
    bundle.createTextDataset("text", text.str());
    return bundlePath;
}

////////////////////////////////////////////////////////////////////

void TextOutFile::addColumn(const string& quantityDescription, string unitDescription, char format, int precision)
{
    _formats.push_back(format);
    _precisions.push_back(precision);

    if (unitDescription.empty()) unitDescription = "1";
    writeLine("# column " + std::to_string(++_ncolumns) + ": " + quantityDescription + " (" + unitDescription + ")");
}

////////////////////////////////////////////////////////////////////

void TextOutFile::writeLine(const string& line)
{
    if (_out.is_open())
    {
        _out << line << "\n";
    }
    else if (!_bundlePath.empty())
    {
        _text += line;
        _text += "\n";
    }
}

////////////////////////////////////////////////////////////////////

void TextOutFile::writeRow(vector<double> values)
{
    writeRowPrivate(values.size(), values.data());
}

////////////////////////////////////////////////////////////////////

void TextOutFile::writeRowPrivate(size_t n, const double* values)
{
    if (n != _ncolumns) throw FATALERROR("Number of values in row does not match the number of columns");

    string line;
    for (size_t i = 0; i < _ncolumns; i++)
    {
        line += (i ? " " : "") + StringUtils::toString(values[i], _formats[i], _precisions[i]);
    }
    writeLine(line);
}

////////////////////////////////////////////////////////////////////
