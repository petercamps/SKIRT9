/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "H5ObjectR.hpp"
#include "H5Core.hpp"

////////////////////////////////////////////////////////////////////

H5ObjectR::H5ObjectR(std::int64_t handle, const char* kind, string path)
    : _handle(handle), _kind(kind), _path(std::move(path))
{}

////////////////////////////////////////////////////////////////////

H5ObjectR::~H5ObjectR()
{
    H5Core::close(_handle);
}

////////////////////////////////////////////////////////////////////

H5ObjectR::H5ObjectR(H5ObjectR&& other) noexcept
    : _handle(other._handle), _kind(other._kind), _path(std::move(other._path))
{
    other._handle = -1;
}

////////////////////////////////////////////////////////////////////

H5ObjectR& H5ObjectR::operator=(H5ObjectR&& other) noexcept
{
    if (this != &other)
    {
        H5Core::close(_handle);
        _handle = other._handle;
        _kind = other._kind;
        _path = std::move(other._path);
        other._handle = -1;
    }
    return *this;
}

////////////////////////////////////////////////////////////////////

string H5ObjectR::context() const
{
    return string(_kind) + " " + _path;
}

////////////////////////////////////////////////////////////////////

bool H5ObjectR::hasAttribute(const string& name) const
{
    return H5Core::hasAttribute(_handle, name);
}

////////////////////////////////////////////////////////////////////

bool H5ObjectR::getBoolAttribute(const string& name) const
{
    return H5Core::boolAttribute(_handle, name, context());
}

////////////////////////////////////////////////////////////////////

int H5ObjectR::getIntAttribute(const string& name) const
{
    return H5Core::intAttribute(_handle, name, context());
}

////////////////////////////////////////////////////////////////////

double H5ObjectR::getDoubleAttribute(const string& name) const
{
    return H5Core::doubleAttribute(_handle, name, context());
}

////////////////////////////////////////////////////////////////////

string H5ObjectR::getStringAttribute(const string& name) const
{
    return H5Core::stringAttribute(_handle, name, context());
}

////////////////////////////////////////////////////////////////////
