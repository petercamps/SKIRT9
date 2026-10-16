/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef H5DATASETW_HPP
#define H5DATASETW_HPP

#include "H5ObjectW.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of the H5DatasetW class represents a dataset in an HDF5 file opened for writing.
    A dataset is created, and its complete contents are written, by one of the create functions of
    the H5BundleW class, which return an instance of this class. The instance then allows setting
    the attributes of the dataset, through the functions of the H5ObjectW base class. See the
    H5Lib class for an overview of the concepts, and the H5ObjectW base class for the rules on
    object lifetime.

    All functions of this class are thread-safe. */
class H5DatasetW final : public H5ObjectW
{
    friend class H5BundleW;

    //================= Construction - Destruction =================

private:
    /** This constructor takes ownership of the specified handle for an open HDF5 dataset. The
        second argument specifies the location of the dataset. */
    H5DatasetW(std::int64_t handle, string path);
};

////////////////////////////////////////////////////////////////////

#endif
