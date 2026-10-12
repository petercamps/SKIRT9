/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef H5SUITER_HPP
#define H5SUITER_HPP

#include "H5BundleR.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of the H5SuiteR class represents a suite in an HDF5 file opened for reading. A
    suite is an HDF5 group that holds a set of bundles, for example the input files for a single
    simulation. Instances of this class are obtained from the H5Lib::openSuite() function. See the
    H5Lib class for an overview of the concepts, and the H5ObjectR base class for the rules on
    object lifetime.

    All functions of this class are thread-safe. */
class H5SuiteR final : public H5ObjectR
{
    friend class H5Lib;

    //================= Construction - Destruction =================

private:
    /** This constructor takes ownership of the specified handle for an open HDF5 group. The
        second argument specifies the location of the group. */
    H5SuiteR(std::int64_t handle, string path);

    //======================== Other Functions =======================

public:
    /** This function returns the names of the bundles located directly in this suite, i.e. the
        names of the HDF5 groups directly inside the group representing the suite, in
        alphabetical order. A group nested inside one of these groups is not listed. */
    vector<string> getBundleNames() const;

    /** This function opens the bundle with the specified name in this suite for reading. It
        throws a fatal error if the suite does not contain a bundle with that name. */
    H5BundleR openBundle(const string& name) const;
};

////////////////////////////////////////////////////////////////////

#endif
