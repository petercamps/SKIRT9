/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "H5SuiteR.hpp"
#include "H5Core.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

H5SuiteR::H5SuiteR(std::int64_t handle, string path) : H5ObjectR(handle, "HDF5 suite", std::move(path)) {}

////////////////////////////////////////////////////////////////////

vector<string> H5SuiteR::getBundleNames() const
{
    return H5Core::groupNames(handle(), context());
}

////////////////////////////////////////////////////////////////////

H5BundleR H5SuiteR::openBundle(const string& name) const
{
    // the location of the suite ends with a colon if the suite is the root group of the file
    string separator = StringUtils::endsWith(path(), ":") ? "" : "/";
    return H5BundleR(H5Core::openGroup(handle(), name, context()), path() + separator + name);
}

////////////////////////////////////////////////////////////////////
