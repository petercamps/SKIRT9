/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef H5BUNDLER_HPP
#define H5BUNDLER_HPP

#include "H5DatasetR.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of the H5BundleR class represents a bundle in an HDF5 file opened for reading. A
    bundle is the counterpart of a regular input file: an HDF5 group containing one or more
    datasets with names fixed by the SKIRT data model for the type of file. Instances of this
    class are obtained from the H5Lib::openBundle() or H5SuiteR::openBundle() functions. See the
    H5Lib class for an overview of the concepts, and the H5ObjectR base class for the rules on
    object lifetime and for the functions to query the attributes of the bundle.

    All functions of this class are thread-safe. */
class H5BundleR final : public H5ObjectR
{
    friend class H5Lib;
    friend class H5SuiteR;

    //================= Construction - Destruction =================

private:
    /** This constructor takes ownership of the specified handle for an open HDF5 group. The
        second argument specifies the location of the group. */
    H5BundleR(std::int64_t handle, string path);

    //======================== Other Functions =======================

public:
    /** This function returns true if the bundle contains a dataset with the specified name, and
        false otherwise. */
    bool hasDataset(const string& name) const;

    /** This function opens the dataset with the specified name in this bundle for reading. It
        throws a fatal error if the bundle does not contain a dataset with that name. */
    H5DatasetR getDataset(const string& name) const;
};

////////////////////////////////////////////////////////////////////

#endif
