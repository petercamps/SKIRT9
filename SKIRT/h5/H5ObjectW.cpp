/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "H5ObjectW.hpp"
#include "H5Core.hpp"

////////////////////////////////////////////////////////////////////

H5ObjectW::H5ObjectW(std::int64_t handle, const char* kind, string path)
    : _handle(handle), _kind(kind), _path(std::move(path))
{}

////////////////////////////////////////////////////////////////////

H5ObjectW::~H5ObjectW()
{
    H5Core::close(_handle);
}

////////////////////////////////////////////////////////////////////

H5ObjectW::H5ObjectW(H5ObjectW&& other) noexcept
    : _handle(other._handle), _kind(other._kind), _path(std::move(other._path))
{
    other._handle = -1;
}

////////////////////////////////////////////////////////////////////

H5ObjectW& H5ObjectW::operator=(H5ObjectW&& other) noexcept
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

string H5ObjectW::context() const
{
    return string(_kind) + " " + _path;
}

////////////////////////////////////////////////////////////////////

void H5ObjectW::setBoolAttribute(const string& name, bool value) const
{
    H5Core::setBoolAttribute(_handle, name, value, context());
}

////////////////////////////////////////////////////////////////////

void H5ObjectW::setIntAttribute(const string& name, int value) const
{
    H5Core::setIntAttribute(_handle, name, value, context());
}

////////////////////////////////////////////////////////////////////

void H5ObjectW::setDoubleAttribute(const string& name, double value) const
{
    H5Core::setDoubleAttribute(_handle, name, value, context());
}

////////////////////////////////////////////////////////////////////

void H5ObjectW::setStringAttribute(const string& name, const string& value) const
{
    H5Core::setStringAttribute(_handle, name, value, context());
}

////////////////////////////////////////////////////////////////////
