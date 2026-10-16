/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef H5LIB_HPP
#define H5LIB_HPP

#include "H5BundleW.hpp"
#include "H5SuiteR.hpp"

////////////////////////////////////////////////////////////////////

/** The H5Lib class offers the entry points for accessing HDF5 files through a thin C++ wrapper
    around the HDF5 C library, offering just the features used by the code. All code outside of
    the \c h5 build target depends only on this wrapper, never directly on the HDF5 library.

    <b>Build option</b>

    Support for HDF5 is optional. It is included when the \c BUILD_HDF5 build option is enabled,
    which requires the HDF5 C library to be installed on the system. Otherwise the code is built
    without any dependency on HDF5. In that case, the available() function returns false, and the
    functions that open or create a suite or a bundle throw a fatal error.

    <b>Concepts</b>

    An HDF5 file is organized like a small file system: it has a root group, which can contain
    further groups as well as datasets, the actual blocks of stored data. Both groups and datasets
    can carry attributes, small named pieces of metadata. The wrapper represents the following
    concepts, each by its own class:

    - A \em bundle (H5BundleR) is the counterpart of a regular input file. It is an HDF5 group
      containing one or more datasets, with names fixed by the SKIRT data model for the type of
      file. For example, a bundle that corresponds to a text column file contains a dataset for
      each column.

    - A \em dataset (H5DatasetR) is a rectangular, multidimensional array of values of the same
      type, located in a bundle.

    - A \em suite (H5SuiteR) is an HDF5 group that holds a set of bundles, for example the input
      files for a single simulation. This allows storing the input for several simulations in a
      single HDF5 file. A suite can be the root group of the file, or a group at any nesting level.

    The names of these classes are suffixed with "R" for reading. A bundle and a dataset also have
    a class for writing, suffixed with "W" (H5BundleW and H5DatasetW, with the common base class
    H5ObjectW). A suite has no class for writing, because writing always targets a specific bundle,
    never a suite as a whole; the suites leading up to a bundle are created as needed.

    <b>Writing</b>

    The createBundle() function opens an HDF5 file for writing, creating the file if it does not
    exist, and creates a new, empty bundle in it, replacing any existing bundle with the same path;
    the rest of the file remains untouched. The datasets of the bundle are created and written
    through the H5BundleW class. A new bundle automatically receives three standard attributes that
    describe its origin: \c producer, the string specified through the setProducer() function;
    \c format, the version of the bundle format (currently 1); and \c created, the date and time at
    which the bundle was created, in ISO 8601 format with milliseconds.

    Only a single process may write to an HDF5 file at any given time, and while it does, no other
    process may read from it. A process may read from and write to the same HDF5 file, but not at
    the same time: all objects obtained for reading must have been destroyed before the file is
    opened for writing, and vice versa.

    <b>Paths</b>

    None of the classes represents an open file; managing the open files is an internal concern
    of the wrapper. Instead, the functions of this class take a path that identifies both the HDF5
    file and the suite or bundle inside it, in the form

        <filepath>:<suite>/<bundle>

    where <tt>\<filepath\></tt> is the path of the HDF5 file, which must have the <tt>.hdf5</tt>
    filename extension, and the portion after the colon is the path of the group inside the file,
    with nested groups separated by forward slashes. The suite portion can be empty, in which case
    the bundle is located in the root group of the file, or it can consist of several nested
    groups. For example:

        /data/campaign7.hdf5:galaxy042/stars.txt

    The path is split at the first occurrence of the <tt>.hdf5</tt> filename extension that is
    followed by a colon. The name of a suite or bundle can therefore end with that extension as
    well. A path without such an occurrence must end with the extension; it refers to the root
    group of the file.

    Combining the name of an input file with the input path and the HDF5 file specified by the
    user to construct such a path, and deciding whether to use a regular file or a bundle, is
    outside the scope of the wrapper. This is handled by the FilePaths class and its clients.

    <b>Object lifetime and thread safety</b>

    Each object handed out by the wrapper acquires its underlying HDF5 resource on construction
    and releases it in its destructor. The objects can be moved but not copied. An object remains
    fully valid after the object from which it was obtained has been destroyed. See the H5ObjectR
    class for more information.

    The HDF5 library cannot be assumed to be thread-safe on its own, since that depends on how
    the installed library was built. All calls into the library are therefore protected by a
    single mutex internal to the wrapper, so that all functions of the wrapper are thread-safe.
    Any number of independent processes may open and read the same HDF5 file at the same time. */
class H5Lib final
{
public:
    /** This function returns true if the code was built with HDF5 support, and false otherwise. */
    static bool available();

    /** This function opens the suite at the specified path for reading. The path has the form
        <tt>\<filepath\>:\<suite\></tt> as described in the class header; if the suite portion is
        empty (with or without the colon), the function opens the root group of the file. The
        function throws a fatal error if the code was built without HDF5 support, if the path is
        malformed, or if the file or the suite cannot be opened. */
    static H5SuiteR openSuite(const string& path);

    /** This function opens the bundle at the specified path for reading, without going through
        its suite first. The path has the form <tt>\<filepath\>:\<suite\>/\<bundle\></tt> as
        described in the class header. The function throws a fatal error if the code was built
        without HDF5 support, if the path is malformed, or if the file or the bundle cannot be
        opened. */
    static H5BundleR openBundle(const string& path);

    /** This function returns true if there is a bundle at the specified path, i.e. if the HDF5
        file exists and holds a group at the location of the bundle, and false if the file does
        not exist or does not hold such a group. The path has the same form as for the
        openBundle() function. The function throws a fatal error if the code was built without
        HDF5 support, if the path is malformed, or if the file exists but cannot be opened as an
        HDF5 file. It allows a caller to report a missing input file in its own words, for example
        to mention a regular file that is missing as well. */
    static bool hasBundle(const string& path);

    /** This function creates a new, empty bundle at the specified path for writing, and returns
        the corresponding H5BundleW instance. The path has the same form as for the openBundle()
        function. If the HDF5 file does not exist, it is created; otherwise it is opened, and an
        existing bundle with the same path is removed, leaving the rest of the file untouched. Any
        missing suites leading up to the bundle are created as well. Note that the HDF5 library
        does not release the space occupied by a removed bundle.

        The new bundle receives the standard attributes \c producer, \c format, and \c created
        described in the class header.

        The function throws a fatal error if the code was built without HDF5 support, if the path
        is malformed, if the file cannot be opened or created for writing (for example because it
        is not an HDF5 file, or because it is in use by another process), if the path refers to a
        dataset or to a suite (a group that contains other groups), or if one of the suites leading
        up to the bundle is a dataset or a bundle. */
    static H5BundleW createBundle(const string& path);

    /** This function sets the value of the \c producer attribute given to every bundle created
        afterwards, normally the name and version of the program. The default value is "SKIRT". The
        function should be called once, early in the program, before any bundle is created. */
    static void setProducer(const string& producer);
};

////////////////////////////////////////////////////////////////////

#endif
