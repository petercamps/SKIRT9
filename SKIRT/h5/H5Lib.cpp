/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "H5Lib.hpp"
#include "FatalError.hpp"
#include "H5Core.hpp"
#include "StringUtils.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // the filename extension that identifies an HDF5 file in a path
    constexpr char EXTENSION[] = ".hdf5";

    // the version of the bundle format, stored as the format attribute of every new bundle
    constexpr int FORMAT = 1;

    // the producer string stored as the producer attribute of every new bundle
    string _producer = "SKIRT";

    // splits the specified path of the form <filepath>:<objectpath> into its two portions, removing any leading and
    // trailing slashes from the object path; the path is split at the first occurrence of the filename extension
    // followed by a colon, so that the name of a group inside the file can end with that extension as well; without
    // such an occurrence, the path must end with the filename extension, and the object path is empty; throws a
    // fatal error if the path is malformed
    std::pair<string, string> splitPath(const string& path)
    {
        string filePath = path;
        string objectPath;
        if (auto pos = path.find(string(EXTENSION) + ":"); pos != string::npos)
        {
            pos += string(EXTENSION).length();
            filePath = path.substr(0, pos);
            objectPath = path.substr(pos + 1);
        }
        else if (!StringUtils::endsWith(path, EXTENSION))
        {
            throw FATALERROR("Path does not refer to an HDF5 file: " + path);
        }

        while (StringUtils::startsWith(objectPath, "/")) objectPath.erase(0, 1);
        while (StringUtils::endsWith(objectPath, "/")) objectPath.pop_back();
        return {filePath, objectPath};
    }
}

////////////////////////////////////////////////////////////////////

bool H5Lib::available()
{
    return H5Core::available();
}

////////////////////////////////////////////////////////////////////

H5SuiteR H5Lib::openSuite(const string& path)
{
    auto [filePath, suitePath] = splitPath(path);
    return H5SuiteR(H5Core::openGroup(filePath, suitePath), filePath + ":" + suitePath);
}

////////////////////////////////////////////////////////////////////

H5BundleR H5Lib::openBundle(const string& path)
{
    auto [filePath, bundlePath] = splitPath(path);
    if (bundlePath.empty()) throw FATALERROR("Path does not refer to a bundle in an HDF5 file: " + path);
    return H5BundleR(H5Core::openGroup(filePath, bundlePath), filePath + ":" + bundlePath);
}

////////////////////////////////////////////////////////////////////

bool H5Lib::hasBundle(const string& path)
{
    auto [filePath, bundlePath] = splitPath(path);
    if (bundlePath.empty()) throw FATALERROR("Path does not refer to a bundle in an HDF5 file: " + path);
    return H5Core::hasGroup(filePath, bundlePath);
}

////////////////////////////////////////////////////////////////////

H5BundleW H5Lib::createBundle(const string& path)
{
    auto [filePath, bundlePath] = splitPath(path);
    if (bundlePath.empty()) throw FATALERROR("Path does not refer to a bundle in an HDF5 file: " + path);

    H5BundleW bundle(H5Core::createGroup(filePath, bundlePath), filePath + ":" + bundlePath);
    bundle.setStringAttribute("producer", _producer);
    bundle.setIntAttribute("format", FORMAT);
    bundle.setStringAttribute("created", System::timestamp(true));
    return bundle;
}

////////////////////////////////////////////////////////////////////

void H5Lib::setProducer(const string& producer)
{
    _producer = producer;
}

////////////////////////////////////////////////////////////////////
