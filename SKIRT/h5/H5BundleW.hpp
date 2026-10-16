/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef H5BUNDLEW_HPP
#define H5BUNDLEW_HPP

#include "Array.hpp"
#include "H5DatasetW.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of the H5BundleW class represents a bundle in an HDF5 file opened for writing. A
    bundle is the counterpart of a regular output file: an HDF5 group containing one or more
    datasets with names fixed by the SKIRT data model for the type of file. Instances of this
    class are obtained from the H5Lib::createBundle() function, which creates a new, empty bundle.
    See the H5Lib class for an overview of the concepts, and the H5ObjectW base class for the rules
    on object lifetime and for the functions to set the attributes of the bundle.

    A dataset is created and written in a single operation, by one of the create functions of this
    class. There is a function for each supported content type: Boolean, 32-bit integer, 64-bit
    floating point, and text. The data is stored in that type, without any conversion. Each of
    these functions returns an H5DatasetW instance, which allows setting the attributes of the new
    dataset.

    The numeric create functions accept the values as a flat sequence and an optional shape, i.e.
    the number of elements along each of the dimensions of the dataset, with the index along the
    last dimension varying most rapidly in the sequence. If the shape is omitted, the dataset is
    one-dimensional, with the number of elements in the sequence. Otherwise the product of the
    dimensions must equal that number. A dimension can be zero, for a dataset without values.

    All functions of this class are thread-safe. */
class H5BundleW final : public H5ObjectW
{
    friend class H5Lib;

    //================= Construction - Destruction =================

private:
    /** This constructor takes ownership of the specified handle for an open HDF5 group. The
        second argument specifies the location of the group. */
    H5BundleW(std::int64_t handle, string path);

    //======================== Other Functions =======================

public:
    /** This function creates a dataset with the specified name and shape in this bundle, and
        writes the specified Boolean values to it. Each value is a character holding 1 for true or
        0 for false; any other nonzero value also means true. The values are stored in the
        representation of a Boolean used by the \c h5py Python package. The function throws a
        fatal error if the number of values does not match the shape, or if the dataset cannot be
        created, for example because the bundle already holds a dataset with that name. */
    H5DatasetW createBoolDataset(const string& name, const vector<char>& values,
                                 const vector<size_t>& shape = vector<size_t>()) const;

    /** This function creates a dataset with the specified name and shape in this bundle, and
        writes the specified values to it as 32-bit signed integers, as described for the
        createBoolDataset() function. */
    H5DatasetW createIntDataset(const string& name, const vector<int>& values,
                                const vector<size_t>& shape = vector<size_t>()) const;

    /** This function creates a dataset with the specified name and shape in this bundle, and
        writes the specified values to it as 64-bit floating point numbers, as described for the
        createBoolDataset() function. */
    H5DatasetW createDoubleDataset(const string& name, const Array& values,
                                   const vector<size_t>& shape = vector<size_t>()) const;

    /** This function creates a dataset with the specified name in this bundle, holding the
        specified UTF-8 encoded text as a single string of variable length, for example the
        complete contents of a text file. The function throws a fatal error if the dataset cannot
        be created. */
    H5DatasetW createTextDataset(const string& name, const string& text) const;

private:
    /** This function returns the shape for a dataset with the specified number of values and the
        shape specified by the caller, which may be empty. It throws a fatal error if the number of
        values does not match the specified shape. The first argument is the name of the dataset,
        for use in error messages. */
    vector<size_t> verifiedShape(const string& name, size_t numValues, const vector<size_t>& shape) const;
};

////////////////////////////////////////////////////////////////////

#endif
