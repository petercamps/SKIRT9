/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef H5OBJECTW_HPP
#define H5OBJECTW_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** H5ObjectW is the common base class for the classes that represent an object in an HDF5 file
    opened for writing: a bundle (H5BundleW) or a dataset (H5DatasetW). See the H5Lib class for an
    overview of these concepts.

    An instance of this class acquires its underlying HDF5 resource on construction and releases
    it in its destructor. Instances can be moved but not copied. An instance remains fully valid
    after the object from which it was obtained has been destroyed. The HDF5 library keeps a file
    open for as long as there is at least one object referring to it. The client is thus
    responsible for eventually destroying every object it obtained; until then, the file remains
    open, and the data written to it may not yet be complete on disk.

    The class offers functions to set the attributes of the object. An attribute is a small, named
    piece of metadata attached to the object. The supported attribute types are Boolean, 32-bit
    integer, 64-bit floating point, and string, each with its own function. Setting an attribute
    that already exists replaces it. The value is stored in the type of the function called,
    without any conversion: a Boolean as an enumeration with the members \c FALSE and \c TRUE (the
    representation used by the \c h5py Python package), an integer as a 32-bit signed integer, a
    floating point number as a 64-bit floating point number, and a string as a string of variable
    length in UTF-8 encoding.

    All functions of this class are thread-safe. */
class H5ObjectW
{
    //================= Construction - Destruction =================

protected:
    /** This constructor takes ownership of the specified handle for an open HDF5 object. The
        second argument names the kind of object (bundle or dataset), and the third argument
        specifies the location of the object, in the form <tt>filepath:objectpath</tt>. Both are
        used in error messages. */
    H5ObjectW(std::int64_t handle, const char* kind, string path);

    /** The destructor releases the underlying HDF5 resource, if any. */
    ~H5ObjectW();

public:
    /** The copy constructor is deleted because instances of this class cannot be copied. */
    H5ObjectW(const H5ObjectW&) = delete;

    /** The copy assignment operator is deleted because instances of this class cannot be copied. */
    H5ObjectW& operator=(const H5ObjectW&) = delete;

    /** The move constructor transfers ownership of the underlying HDF5 resource from the
        specified instance to the newly constructed one. The specified instance is left in a state
        in which it can only be destroyed or assigned to. */
    H5ObjectW(H5ObjectW&& other) noexcept;

    /** The move assignment operator releases the underlying HDF5 resource of the receiving
        instance, if any, and transfers ownership of the resource of the specified instance to the
        receiving instance. The specified instance is left in a state in which it can only be
        destroyed or assigned to. */
    H5ObjectW& operator=(H5ObjectW&& other) noexcept;

    //======================== Other Functions =======================

public:
    /** This function sets the Boolean attribute with the specified name to the specified value.
        */
    void setBoolAttribute(const string& name, bool value) const;

    /** This function sets the integer attribute with the specified name to the specified value,
        stored as a 32-bit signed integer. */
    void setIntAttribute(const string& name, int value) const;

    /** This function sets the floating point attribute with the specified name to the specified
        value, stored as a 64-bit floating point number. */
    void setDoubleAttribute(const string& name, double value) const;

    /** This function sets the string attribute with the specified name to the specified UTF-8
        encoded value. */
    void setStringAttribute(const string& name, const string& value) const;

protected:
    /** This function returns the handle for the underlying HDF5 object, for use by subclasses. */
    std::int64_t handle() const { return _handle; }

    /** This function returns the location of the object that was specified at construction, for
        use by subclasses. */
    const string& path() const { return _path; }

    /** This function returns a human-readable description of the object for use in error
        messages, including the kind and the location of the object. */
    string context() const;

    //======================== Data Members ========================

private:
    std::int64_t _handle{-1};
    const char* _kind{""};
    string _path;
};

////////////////////////////////////////////////////////////////////

#endif
