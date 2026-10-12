/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef H5DATASETR_HPP
#define H5DATASETR_HPP

#include "Array.hpp"
#include "H5ObjectR.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of the H5DatasetR class represents a dataset in an HDF5 file opened for reading,
    i.e. a rectangular, multidimensional array of values of the same type. Instances of this class
    are obtained from the H5BundleR::getDataset() function. See the H5Lib class for an overview of
    the concepts, and the H5ObjectR base class for the rules on object lifetime and for the
    functions to query the attributes of the dataset.

    The supported content types are Boolean, 32-bit integer, and 64-bit floating point, each with
    its own read function. The contents are returned in the type of the function called,
    regardless of how they are stored in the file. The HDF5 library performs the conversion (e.g.
    between 32-bit and 64-bit number representations) and transparently handles chunked or
    compressed data. If a stored value cannot be represented in that type without loss of
    information, a fatal error is thrown.

    All functions of this class are thread-safe. */
class H5DatasetR final : public H5ObjectR
{
    friend class H5BundleR;

    //================= Construction - Destruction =================

private:
    /** This constructor takes ownership of the specified handle for an open HDF5 dataset. The
        second argument specifies the location of the dataset. */
    H5DatasetR(std::int64_t handle, string path);

    //======================== Other Functions =======================

public:
    /** This function returns the shape of the dataset, i.e. the number of elements along each of
        its dimensions. The dimensions are listed in the order in which they were specified when
        the dataset was created; the index along the last dimension varies most rapidly in the
        sequence of values returned by the read functions. The result is empty for a dataset
        holding a single scalar value. */
    vector<size_t> shape() const;

    /** This function returns the complete contents of the dataset as Boolean values. Each element
        of the result is 1 for true or 0 for false; the result is a vector of \c char rather than
        \c bool because the latter is a specialization that does not behave like a regular
        container. The values are returned as a flat sequence with the index along the last
        dimension of the dataset varying most rapidly; use the shape() function to obtain the
        dimensions. The values must be stored as an integer or enumeration type; a nonzero stored
        value means true. Otherwise the function throws a fatal error. */
    vector<char> readBool() const;

    /** This function returns the complete contents of the dataset as 32-bit signed integers. The
        values are returned as a flat sequence with the index along the last dimension of the
        dataset varying most rapidly; use the shape() function to obtain the dimensions. The
        values must be stored as an integer or floating point type, and each value must be
        representable as a 32-bit signed integer without loss of information. Otherwise the
        function throws a fatal error. */
    vector<int> readInt() const;

    /** This function returns the complete contents of the dataset as 64-bit floating point
        numbers. The values are returned as a flat sequence with the index along the last
        dimension of the dataset varying most rapidly; use the shape() function to obtain the
        dimensions. The values must be stored as an integer or floating point type, and each value
        must be representable as a 64-bit floating point number without loss of information.
        Otherwise the function throws a fatal error. */
    Array readDouble() const;
};

////////////////////////////////////////////////////////////////////

#endif
