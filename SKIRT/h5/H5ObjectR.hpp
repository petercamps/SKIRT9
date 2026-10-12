/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef H5OBJECTR_HPP
#define H5OBJECTR_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** H5ObjectR is the common base class for the classes that represent an object in an HDF5 file
    opened for reading: a suite (H5SuiteR), a bundle (H5BundleR), or a dataset (H5DatasetR). See
    the H5Lib class for an overview of these concepts.

    An instance of this class acquires its underlying HDF5 resource on construction and releases
    it in its destructor. Instances can be moved but not copied. An instance remains fully valid
    after the object from which it was obtained has been destroyed; for example, a dataset does not
    depend on its bundle still being alive. The HDF5 library keeps a file open for as long as there
    is at least one object referring to it. The client is thus responsible for eventually
    destroying every object it obtained; until then, the file remains open.

    The class offers functions to query the attributes of the object. An attribute is a small,
    named piece of metadata attached to the object. The supported attribute types are Boolean,
    32-bit integer, 64-bit floating point, and string, each with its own function. The value is
    returned in the type of the function called, regardless of how it is stored in the file, as
    long as the stored value can be represented in that type without loss of information;
    otherwise a fatal error is thrown.

    All functions of this class are thread-safe. */
class H5ObjectR
{
    //================= Construction - Destruction =================

protected:
    /** This constructor takes ownership of the specified handle for an open HDF5 object. The
        second argument names the kind of object (suite, bundle, or dataset), and the third
        argument specifies the location of the object, in the form <tt>filepath:objectpath</tt>.
        Both are used in error messages. */
    H5ObjectR(std::int64_t handle, const char* kind, string path);

    /** The destructor releases the underlying HDF5 resource, if any. */
    ~H5ObjectR();

public:
    /** The copy constructor is deleted because instances of this class cannot be copied. */
    H5ObjectR(const H5ObjectR&) = delete;

    /** The copy assignment operator is deleted because instances of this class cannot be copied. */
    H5ObjectR& operator=(const H5ObjectR&) = delete;

    /** The move constructor transfers ownership of the underlying HDF5 resource from the
        specified instance to the newly constructed one. The specified instance is left in a state
        in which it can only be destroyed or assigned to. */
    H5ObjectR(H5ObjectR&& other) noexcept;

    /** The move assignment operator releases the underlying HDF5 resource of the receiving
        instance, if any, and transfers ownership of the resource of the specified instance to the
        receiving instance. The specified instance is left in a state in which it can only be
        destroyed or assigned to. */
    H5ObjectR& operator=(H5ObjectR&& other) noexcept;

    //======================== Other Functions =======================

public:
    /** This function returns true if the object has an attribute with the specified name, and
        false otherwise. It is needed because in some places the absence of an attribute is
        meaningful, and distinct from an empty value. */
    bool hasAttribute(const string& name) const;

    /** This function returns the value of the Boolean attribute with the specified name. The
        attribute must hold a single value stored as an integer or enumeration type; a nonzero
        value means true. The function throws a fatal error if the attribute does not exist or if
        its value cannot be returned as a Boolean. */
    bool getBoolAttribute(const string& name) const;

    /** This function returns the value of the integer attribute with the specified name. The
        attribute must hold a single value stored as an integer or floating point type, which can
        be represented as a 32-bit signed integer without loss of information. The function throws
        a fatal error if the attribute does not exist or if its value cannot be returned as an
        integer. */
    int getIntAttribute(const string& name) const;

    /** This function returns the value of the floating point attribute with the specified name.
        The attribute must hold a single value stored as an integer or floating point type. The
        function throws a fatal error if the attribute does not exist or if its value cannot be
        returned as a floating point number. */
    double getDoubleAttribute(const string& name) const;

    /** This function returns the value of the string attribute with the specified name, as a
        UTF-8 encoded string. The attribute must hold a single string of fixed or variable length.
        The function throws a fatal error if the attribute does not exist or if its value is not a
        string. */
    string getStringAttribute(const string& name) const;

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
