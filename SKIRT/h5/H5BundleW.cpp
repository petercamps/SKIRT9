/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "H5BundleW.hpp"
#include "FatalError.hpp"
#include "H5Core.hpp"

////////////////////////////////////////////////////////////////////

H5BundleW::H5BundleW(std::int64_t handle, string path) : H5ObjectW(handle, "HDF5 bundle", std::move(path)) {}

////////////////////////////////////////////////////////////////////

vector<size_t> H5BundleW::verifiedShape(const string& name, size_t numValues, const vector<size_t>& shape) const
{
    if (shape.empty()) return vector<size_t>{numValues};

    size_t product = 1;
    for (size_t size : shape) product *= size;
    if (product != numValues)
        throw FATALERROR("The number of values for dataset '" + name + "' in " + context()
                         + " does not match its shape");
    return shape;
}

////////////////////////////////////////////////////////////////////

H5DatasetW H5BundleW::createBoolDataset(const string& name, const vector<char>& values,
                                        const vector<size_t>& shape) const
{
    // reduce each value to one for true or zero for false
    vector<char> flags(values.size());
    for (size_t i = 0; i != values.size(); ++i) flags[i] = values[i] != 0;

    auto dims = verifiedShape(name, flags.size(), shape);
    return H5DatasetW(H5Core::createBoolDataset(handle(), name, dims, flags.data(), context()), path() + "/" + name);
}

////////////////////////////////////////////////////////////////////

H5DatasetW H5BundleW::createIntDataset(const string& name, const vector<int>& values, const vector<size_t>& shape) const
{
    auto dims = verifiedShape(name, values.size(), shape);
    return H5DatasetW(H5Core::createIntDataset(handle(), name, dims, values.data(), context()), path() + "/" + name);
}

////////////////////////////////////////////////////////////////////

H5DatasetW H5BundleW::createDoubleDataset(const string& name, const Array& values, const vector<size_t>& shape) const
{
    auto dims = verifiedShape(name, values.size(), shape);
    return H5DatasetW(H5Core::createDoubleDataset(handle(), name, dims, begin(values), context()), path() + "/" + name);
}

////////////////////////////////////////////////////////////////////

H5DatasetW H5BundleW::createTextDataset(const string& name, const string& text) const
{
    return H5DatasetW(H5Core::createTextDataset(handle(), name, text, context()), path() + "/" + name);
}

////////////////////////////////////////////////////////////////////
