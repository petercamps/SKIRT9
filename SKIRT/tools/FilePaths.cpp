/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FilePaths.hpp"
#include "FatalError.hpp"
#include "H5Lib.hpp"
#include "ProcessManager.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <mutex>
#include <unordered_map>

////////////////////////////////////////////////////////////////////

namespace
{
    // flag becomes true if the static paths have been initialized
    std::once_flag _initialized;

    // the resource paths:  <filename, complete_path>
    std::unordered_map<string, string> _resourcePaths;

    // the expected resource pack names
    vector<string> _expectedPacks;

    // the expected resource pack version numbers:  <packname, version>
    std::unordered_map<string, int> _expectedPackVersions;

    // the installed resource pack version numbers:  <packname, version>
    std::unordered_map<string, int> _installedPackVersions;

    // relative paths to check for presence of built-in resources
    constexpr const char* _intpaths[] = {"../../../git/SKIRT/resources", "../../../../git/SKIRT/resources"};
    const int _Nintpaths = sizeof(_intpaths) / sizeof(const char*);

    // relative paths to check for presence of external resources
    constexpr const char* _extpaths[] = {"../../../resources", "../../../../resources"};
    const int _Nextpaths = sizeof(_extpaths) / sizeof(const char*);

    // recursively searches the given directory and its subdirectories for any files,
    // adds the corresponding canonical paths to the resource dictionary,
    // and processes any "version.txt" files encountered
    void findResourcesIn(const string& directory)
    {
        // search the current level
        for (const string& filename : System::filesInDirectory(directory))
        {
            if (!_resourcePaths.count(filename))
                _resourcePaths.emplace(filename, System::canonicalPath(StringUtils::joinPaths(directory, filename)));

            // remember the name and version number for each installed resource pack
            if (filename == "version.txt")
            {
                auto segments = StringUtils::split(directory, "_");
                if (segments.size() > 1)
                {
                    string packname = segments[segments.size() - 1];
                    int packversion = 0;
                    try
                    {
                        std::ifstream versionfile(StringUtils::joinPaths(directory, filename));
                        versionfile >> packversion;
                    }
                    catch (...)
                    {};

                    if (packversion > 0 && !_installedPackVersions.count(packname))
                        _installedPackVersions.emplace(packname, packversion);
                }
            }
        }

        // search the subdirectories
        for (const string& subdirname : System::dirsInDirectory(directory))
        {
            findResourcesIn(StringUtils::joinPaths(directory, subdirname));
        }
    }

    // parses the list of expected resource packs, if present, and stores the results
    // !! assumes that the dictionary with resource paths has already been populated !!
    void readExpectedPacks()
    {
        if (_resourcePaths.count("ExpectedResources.txt"))
        {
            try
            {
                std::ifstream expectedfile(_resourcePaths.at("ExpectedResources.txt"));
                while (expectedfile.good())
                {
                    string packname;
                    int packversion = 0;
                    expectedfile >> packname >> packversion;
                    if (!packname.empty() && packversion > 0 && !_expectedPackVersions.count(packname))
                    {
                        _expectedPacks.emplace_back(packname);
                        _expectedPackVersions.emplace(packname, packversion);
                    }
                }
            }
            catch (...)
            {};
        }
    }

    // populates the dictionary with resource paths and gathers resource pack information,
    // or throws a fatal error if there is a problem
    void findResources()
    {
        // get the executable path (or the empty string in case of failure)
        string executableFilePath = System::executablePath();
        if (executableFilePath.empty()) throw FATALERROR("Could not determine path to executable");

        // get the location of the executable (i.e. the path to the containing directory)
        string executableDirPath = StringUtils::dirPath(executableFilePath);

        // iterate over the relative paths for built-in resources
        for (int i = 0; i < _Nintpaths; i++)
        {
            string test = StringUtils::joinPaths(executableDirPath, _intpaths[i]);
            if (System::isDir(test)) findResourcesIn(test);
        }

        // at this point we should have found at least some internal resources
        if (_resourcePaths.empty())
            throw FATALERROR("Could not locate built-in resources relative to '" + executableDirPath + "'");

        // iterate over the relative paths for external resources
        for (int i = 0; i < _Nextpaths; i++)
        {
            string test = StringUtils::joinPaths(executableDirPath, _extpaths[i]);
            if (System::isDir(test)) findResourcesIn(test);
        }

        // parse the list of expected resource packs, if present
        readExpectedPacks();
    }

    // returns true if the resource filename matches the requirements, false otherwise;
    // see documentation of the resourceName() function for more details
    bool matches(const string& resource, const string& type, const vector<string>& segments)
    {
        if (!StringUtils::endsWith(resource, type)) return false;
        for (const string& segment : segments)
        {
            if (!StringUtils::contains(resource, segment)) return false;
        }
        return true;
    }

    // returns a human-readable message describing the specified requirements;
    // see documentation of the resourceName() function for more details
    string message(const string& type, const vector<string>& segments)
    {
        string msg = "type '" + type + "'";
        if (!segments.empty())
        {
            msg += " with filename containing";
            for (const string& segment : segments) msg += " '" + segment + "',";
            msg.erase(msg.size() - 1, 1);
        }
        return msg;
    }
}

////////////////////////////////////////////////////////////////////

FilePaths::FilePaths(SimulationItem* parent)
{
    parent->addChild(this);
}

////////////////////////////////////////////////////////////////////

void FilePaths::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();
    std::call_once(_initialized, findResources);
}

////////////////////////////////////////////////////////////////////

namespace
{
    // If the specified path refers to an HDF5 file, this function returns the portion of the path up to and including
    // the name of the HDF5 file, and the portion that follows the colon after that name, with leading and trailing
    // slashes removed. Otherwise, it returns two empty strings. A path refers to an HDF5 file if it contains the HDF5
    // filename extension followed by a colon or, if the second argument is true, if it ends with that extension.
    std::pair<string, string> splitHdf5Reference(const string& path, bool allowAtEnd)
    {
        const string extension = ".hdf5";
        string filePath;
        string objectPath;
        if (auto pos = path.find(extension + ":"); pos != string::npos)
        {
            filePath = path.substr(0, pos + extension.length());
            objectPath = path.substr(pos + extension.length() + 1);
            while (StringUtils::startsWith(objectPath, "/")) objectPath.erase(0, 1);
            while (StringUtils::endsWith(objectPath, "/")) objectPath.pop_back();
        }
        else if (allowAtEnd && StringUtils::endsWith(path, extension))
        {
            filePath = path;
        }
        return {filePath, objectPath};
    }

    // returns the absolute canonical path of the specified directory, ending with a single (forward or backward) slash
    string canonicalDirPath(const string& dirPath)
    {
        string result = System::canonicalPath(dirPath);
        if (!StringUtils::endsWith(result, "/") && !StringUtils::endsWith(result, "\\")) result += "/";
        return result;
    }

    // throws a fatal error if the code was built without HDF5 support; the arguments describe and specify the path
    // that refers to an HDF5 file
    void requireHdf5(const string& description, const string& path)
    {
        if (!H5Lib::available())
            throw FATALERROR(description + " '" + path + "' refers to an HDF5 file, but this version of the code was "
                             + "built without HDF5\n(enable the BUILD_HDF5 build option and rebuild the code)");
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // splits the specified value of the input or output path, of the form <dir>/<hdf>:<suite>, into the path of the
    // directory, the path of the HDF5 file (or the empty string), and the path of the suite (or the empty string),
    // and verifies that the directory exists; the first argument is "Input" or "Output", for use in error messages
    std::tuple<string, string, string> splitLocation(const string& kind, const string& value)
    {
        // unless the value is an existing directory, it may include an HDF5 file and a suite after the directory
        if (System::isDir(value)) return {value, string(), string()};

        auto [hdf5Path, suite] = splitHdf5Reference(value, true);
        if (hdf5Path.empty()) throw FATALERROR(kind + " path does not exist or is not a directory: " + value);
        requireHdf5(kind + " path", value);

        string dirPath = StringUtils::dirPath(hdf5Path);
        if (dirPath.empty()) dirPath = ".";
        if (!System::isDir(dirPath)) throw FATALERROR(kind + " path does not exist or is not a directory: " + dirPath);
        return {dirPath, hdf5Path, suite};
    }
}

////////////////////////////////////////////////////////////////////

void FilePaths::setInputPath(const string& value)
{
    auto [dirPath, hdf5Path, suite] = splitLocation("Input", value);
    if (!hdf5Path.empty() && !System::isFile(hdf5Path)) throw FATALERROR("Input HDF5 file does not exist: " + hdf5Path);

    _inputPath = canonicalDirPath(dirPath);
    _inputHdf5Path = hdf5Path.empty() ? string() : _inputPath + StringUtils::filename(hdf5Path);
    _inputSuite = suite;
}

////////////////////////////////////////////////////////////////////

std::pair<string, string> FilePaths::input(const string& name) const
{
    // if the name embeds a reference to a bundle in an HDF5 file, there is no regular file candidate
    if (auto [hdf5Path, bundlePath] = splitHdf5Reference(name, false); !hdf5Path.empty())
    {
        requireHdf5("Input file name", name);
        if (!StringUtils::isAbsolutePath(hdf5Path)) hdf5Path = _inputPath + hdf5Path;
        return {string(), hdf5Path + ":" + bundlePath};
    }

    // an absolute path refers to a regular file, and has no bundle candidate
    if (StringUtils::isAbsolutePath(name)) return {name, string()};

    // otherwise, the name is relative to the input directory, and to the HDF5 file and suite in the input path, if any
    string bundlePath;
    if (!_inputHdf5Path.empty()) bundlePath = _inputHdf5Path + ":" + StringUtils::joinPaths(_inputSuite, name);
    return {_inputPath + name, bundlePath};
}

////////////////////////////////////////////////////////////////////

bool FilePaths::isAbsolutePath(const string& path)
{
    string hdf5Path = splitHdf5Reference(path, true).first;
    return StringUtils::isAbsolutePath(hdf5Path.empty() ? path : hdf5Path);
}

////////////////////////////////////////////////////////////////////

void FilePaths::setOutputPath(const string& value)
{
    auto [dirPath, hdf5Path, suite] = splitLocation("Output", value);
    if (!hdf5Path.empty() && System::isDir(hdf5Path)) throw FATALERROR("Output HDF5 file is a directory: " + hdf5Path);

    _outputPath = canonicalDirPath(dirPath);
    _outputHdf5Path = hdf5Path.empty() ? string() : _outputPath + StringUtils::filename(hdf5Path);
    _outputSuite = suite;
}

////////////////////////////////////////////////////////////////////

void FilePaths::setOutputPrefix(string value)
{
    _outputPrefix = value;
}

////////////////////////////////////////////////////////////////////

string FilePaths::outputPrefix() const
{
    return _outputPrefix;
}

////////////////////////////////////////////////////////////////////

std::pair<string, string> FilePaths::output(const string& name) const
{
    string filename = _outputPrefix + "_" + name;
    string bundlePath;
    if (!_outputHdf5Path.empty()) bundlePath = _outputHdf5Path + ":" + StringUtils::joinPaths(_outputSuite, filename);
    return {_outputPath + filename, bundlePath};
}

////////////////////////////////////////////////////////////////////

string FilePaths::directoryPath(const string& path)
{
    if (System::isDir(path)) return path;
    string hdf5Path = splitHdf5Reference(path, true).first;
    return hdf5Path.empty() ? path : StringUtils::dirPath(hdf5Path);
}

////////////////////////////////////////////////////////////////////

void FilePaths::setSharedProcess(bool shared)
{
    _sharedProcess = shared;
}

////////////////////////////////////////////////////////////////////

void FilePaths::cutOver()
{
    if (_cutOver) return;

    // a process that runs a single simulation must have released all objects obtained while reading bundles
    if (!_sharedProcess)
    {
        auto files = H5Lib::openFiles();
        if (!files.empty())
            throw FATALERROR("HDF5 files are still open at the end of setup: " + StringUtils::join(files, ", "));
    }

    // the root process verifies that it can write the HDF5 output file, if any
    if (!_outputHdf5Path.empty() && ProcessManager::isRoot()) H5Lib::verifyWritable(_outputHdf5Path);

    // record the cutover only if the checks succeed
    _cutOver = true;
}

////////////////////////////////////////////////////////////////////

bool FilePaths::hasCutOver() const
{
    return _cutOver;
}

////////////////////////////////////////////////////////////////////

void FilePaths::verifyBundleInput(const string& bundlePath) const
{
    if (_cutOver)
        throw FATALERROR("Bundle " + bundlePath
                         + " cannot be read after the end of setup, because all input must be read during setup");
}

////////////////////////////////////////////////////////////////////

void FilePaths::verifyBundleOutput(const string& bundlePath) const
{
    if (!_cutOver)
        throw FATALERROR("Bundle " + bundlePath
                         + " cannot be written before the end of setup, because output must be written after setup");
}

////////////////////////////////////////////////////////////////////

bool FilePaths::hasResource(string name)
{
    std::call_once(_initialized, findResources);

    return _resourcePaths.count(name) > 0;
}

////////////////////////////////////////////////////////////////////

string FilePaths::resource(const string& name)
{
    std::call_once(_initialized, findResources);

    if (!_resourcePaths.count(name))
    {
        throw FATALERROR("Could not locate resource '" + name
                         + "'"
                           "\nRun ./downloadResources.sh in the ~/SKIRT/git directory"
                           "\nOr download additional resources from www.skirt.ugent.be");
    }
    return _resourcePaths.at(name);
}

////////////////////////////////////////////////////////////////////

string FilePaths::resourceName(const string& type, const vector<string>& segments)
{
    std::call_once(_initialized, findResources);

    // look for matching resource filename
    string result;
    for (const auto& [resource, path] : _resourcePaths)
    {
        if (matches(resource, type, segments))
        {
            // fail if there is ambiguity
            if (!result.empty()) throw FATALERROR("Multiple resources matching " + message(type, segments));
            result = resource;
        }
    }
    // fail if not found
    if (result.empty()) throw FATALERROR("No resources matching " + message(type, segments));
    return result;
}

////////////////////////////////////////////////////////////////////

vector<string> FilePaths::expectedPacks()
{
    std::call_once(_initialized, findResources);

    return _expectedPacks;
}

////////////////////////////////////////////////////////////////////

int FilePaths::expectedPackVersion(const string& name)
{
    std::call_once(_initialized, findResources);

    return _expectedPackVersions.count(name) ? _expectedPackVersions.at(name) : 0;
}

////////////////////////////////////////////////////////////////////

int FilePaths::installedPackVersion(const string& name)
{
    std::call_once(_initialized, findResources);

    return _installedPackVersions.count(name) ? _installedPackVersions.at(name) : 0;
}

////////////////////////////////////////////////////////////////////
