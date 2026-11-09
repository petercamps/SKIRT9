/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FILEPATHS_HPP
#define FILEPATHS_HPP

#include "SimulationItem.hpp"
#include <unordered_map>

////////////////////////////////////////////////////////////////////

/** The FilePaths class manages the paths for the input and output files of a simulation, and for
    the resources included with the code or provided externally.

    <b>Input and output files</b>

    A client program determines the input and output file paths and the prefix for output file
    names, presumably from the command line options. It then stores this information in the
    FilePaths instance held by the simulation, so that the simulation hierarchy can easily access
    it.

    <b>Input from HDF5 files</b>

    Input files can also be provided as \em bundles in an HDF5 file, as described for the H5Lib
    class. To this end, the input path specified by the client program can include up to three
    components:

        <dir>/<hdf>:<suite>

    where <tt>\<dir\></tt> is the input directory; <tt>\<hdf\></tt> is the name of an HDF5 file
    located in that directory, including the mandatory <tt>.hdf5</tt> filename extension; and
    <tt>\<suite\></tt> is the path of a group in that HDF5 file, with nested groups separated by
    forward slashes. Only the directory is required, and a suite can be specified only together
    with an HDF5 file.

    If only a directory is given, every input file is expected to be a regular file in that
    directory. If an HDF5 file is given as well, an input file is first sought as a regular file
    in the directory, and if there is no such file, as a bundle in the HDF5 file, with a name that
    matches the name of the input file. If a suite is given, the bundle is sought in that suite
    rather than in the root group of the HDF5 file. This makes it possible to store the input for
    several simulations in a single HDF5 file, each in its own suite.

    Furthermore, the name of an input file, as configured for example in a ski file, can itself
    embed a reference to a bundle in an HDF5 file, in the form

        <hdf>:<suite>/<bundle>

    In that case, the input file is always read from the specified bundle, regardless of the HDF5
    file and suite that may be included in the input path. A relative path for the HDF5 file is
    interpreted relative to the input directory.

    The input() function implements these rules by returning two candidate locations for an input
    file, a regular file and a bundle. It is up to the caller to try the regular file first, and to
    fall back to the bundle only if that file does not exist.

    <b>Output to HDF5 files</b>

    In the same way, the output path specified by the client program can include a directory, an
    HDF5 file in that directory, and a suite in that file. If an HDF5 file is given, output files
    that support it are written as bundles in that file instead of as regular files in the
    directory, in the root group of the file or in the specified suite. The HDF5 file is created if
    it does not exist, and otherwise extended; a bundle with the same name as one being written is
    replaced. The bundle name is the name of the regular file that would otherwise have been
    written. Output files that are always written as regular files, such as the log file, are
    placed in the directory.

    The output() and outputLocation() functions implement these rules by returning two candidate
    locations for an output file, a regular file and a bundle. A caller that supports writing a
    bundle uses the bundle candidate if it is nonempty, and the regular file otherwise.

    <b>Duplicate output file names</b>

    Output files are named for the simulation items that write them, so that two items with
    similar names may write output files with the same name, in which case the second file
    overwrites the first (or the second bundle replaces the first). To report this situation, the
    output() function, which a class calls once for each output file it writes, registers the file
    name together with the simulation item on whose behalf the file is written. If the name has
    been registered before, the function logs a warning that names both simulation items. The
    names are compared regardless of case, because names that differ only in case refer to the
    same file on a file system that ignores case, such as the default file systems on macOS and
    Windows; the warning then mentions this condition. A warning is logged only for the first
    conflict between a given pair of items (or for a given item, if it writes a file with the same
    name more than once), so that two items that write a series of files with the same names
    cause a single warning. The registration covers the files written since the FilePaths instance
    was constructed, i.e. by a single simulation. The warning cannot prevent the overwrite, and it
    is usually issued after the simulation has run, since most output files are written at the
    end.

    The outputLocation() function returns the same candidate locations without registering the
    name. It is used to locate the log file and the parameter file, whose names cannot conflict
    with those of the files written by simulation items, which always have the \c dat or \c fits
    filename extension; to locate a file that has already been written, for example to store a
    copy as a bundle; and to determine whether the output path includes an HDF5 file.

    <b>I/O cutover</b>

    The I/O cutover is the moment in a simulation that separates all reading of bundles from all
    writing of bundles. It falls at the end of setup, after all processes have finished their
    setup: every input file is read during setup, and every output file is written after it.
    Before the cutover, bundles may be read but not written; after the cutover, bundles may be
    written but not read, nor tested for existence. Regular files are not affected. Together with
    the fact that only the root process writes bundles, this guarantees that the processes of a
    multiprocessing run never access an HDF5 file in a way that violates the concurrency rules of
    HDF5 (see the H5Lib class).

    The simulation calls the cutOver() function at the cutover. The classes that read or write
    bundles call the verifyBundleInput() or verifyBundleOutput() function before each access, so
    that a violation of the rule is reported as a fatal error in every run, also when there is a
    single process. In a process that runs a single simulation, the cutOver() function also
    verifies that the process holds no open HDF5 objects, so that a reader that keeps a bundle
    open is reported right away. This check is skipped when several simulations run in parallel in
    the same process, as specified through the setSharedProcess() function, because the other
    simulations may legitimately hold HDF5 objects at that moment. Finally, if the output path
    includes an HDF5 file, the cutOver() function verifies in the root process that the file can
    be written, so that a file that is not an HDF5 file, that is read-only, or that is held open by
    another program is reported before the simulation runs, rather than when the first bundle is
    written.

    <b>Resource files</b>

    The FilePaths class also offers static functions for locating the resource files used by the
    simulation hierarchy. These files can be provided as part of the build tree in the source code
    repository, or can be installed in a directory next to (and thus outside of) the build tree.
    This mechanism allows to provide small and frequently-used resource files as part of the source
    code repository, while requiring larger resource files to be downloaded seperately.

    Resource files are identified by their filename, without any directory or path information. In
    other words, a given resource file could be located in any of the supported directories (or
    nested subdirectory). On the other hand, this means there can be only a single resource with a
    particular name.

    Specifically, when first invoked, the FilePaths class builds a list of all available resource
    files by iterating over the following directories and all nested subdirectories inside these
    directories, recursively:

    - the \c resources directory inside the SKIRT build tree.

    - the \c resources directory (if any) next to the SKIRT \c git directory (i.e. outside of the
    build tree).

    The top-level directories are iterated in the order listed above. The iteration order for the
    nested subdirectories inside the top-level directories is unspecified. The first occurrence of
    a particular filename is stored; any subsequent occurrences of the same filename are ignored.

    <b>Resource packs</b>

    Finally, the FilePaths class offers static functions that support versioning for resource \em
    packs, which are usually installed in the resource directory outside of the build tree.

    A resource pack is a set of resource files contained in a resource subdirectory, called the
    pack directory, that contains a file named \c version.txt. The first token in this text file's
    contents is interpreted as the integer version number of the resource pack. The name of the
    resource pack is determined by the segment of the pack directory name after the last
    underscore. The resource files in the pack can reside in nested subdirectories.

    Furthermore, a resource file named \c ExpectedResources.txt can be provided as part of the
    build tree. This text file lists the names and version numbers of the resource packs expected
    to be installed, one pack name and corresponding version number per line.

    The FilePaths class offers functions to retrieve information on the expected and installed
    resource packs from these files. This information can be used by the client program to verify
    that the appropriate resource pack versions have been installed.

    */
class FilePaths : public SimulationItem
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a file path object that is hooked up as a child to the specified
        parent in the simulation hierarchy, so that it will automatically be deleted. The setup()
        function is \em not called by this constructor. */
    explicit FilePaths(SimulationItem* parent);

protected:
    /** This function determines and caches the resource file paths that can be returned by this
        class. This avoids repeated searches through the resource directories, and allows reporting
        any problems as early as possible in the program's lifecycle. */
    void setupSelfBefore() override;

    //======== Input and output file paths and names =======

public:
    /** This function sets the location of the input files. The argument is the (absolute or
        relative) path of the input directory, optionally followed by the name of an HDF5 file in
        that directory and by a suite in that file, in the form <tt>\<dir\>/\<hdf\>:\<suite\></tt>
        described in the class header. An empty string (the default value) means the current
        directory. The function throws a fatal error if the directory or the HDF5 file does not
        exist, or if the argument refers to an HDF5 file while the code was built without HDF5
        support. It does not verify that the HDF5 file can be opened or that the suite exists. */
    void setInputPath(const string& value);

    /** This function determines the candidate locations for the input file with the specified
        name, as described in the class header. It returns a pair of strings.

        The first string is the absolute path of a regular file, formed by interpreting the
        specified name relative to the input directory (or the specified name itself if it is an
        absolute path).

        The second string is the path of a bundle in an HDF5 file, in the form
        <tt>\<filepath\>:\<suite\>/\<bundle\></tt> expected by the H5Lib::openBundle() function.
        It is formed by using the specified name as the name of a bundle in the HDF5 file and
        suite included in the input path, if any.

        Either string is empty if there is no corresponding candidate. Specifically, the second
        string is empty if the input path does not include an HDF5 file, or if the specified name
        is an absolute path. If the specified name embeds a reference to a bundle in an HDF5 file,
        in the form <tt>\<hdf\>:\<suite\>/\<bundle\></tt>, the first string is empty and the
        second string is formed from that reference, with a relative HDF5 file path interpreted
        relative to the input directory. In that case, the function throws a fatal error if the
        code was built without HDF5 support.

        The function performs no input/output operations, so a nonempty string does not guarantee
        that the file or bundle exists. The caller is expected to try the regular file first, and
        to fall back to the bundle only if that file does not exist. */
    std::pair<string, string> input(const string& name) const;

    /** This function returns true if the specified input path is absolute, and false if it is
        relative. The path can have the form <tt>\<dir\>/\<hdf\>:\<suite\></tt> described in the
        class header, in which case only the portion up to and including the name of the HDF5 file
        is considered. In contrast, the StringUtils::isAbsolutePath() function considers any path
        containing a colon to be absolute. */
    static bool isAbsolutePath(const string& path);

    /** This function sets the location of the output files. The argument is the (absolute or
        relative) path of the output directory, optionally followed by the name of an HDF5 file in
        that directory and by a suite in that file, in the form <tt>\<dir\>/\<hdf\>:\<suite\></tt>
        described in the class header. An empty string (the default value) means the current
        directory. The function throws a fatal error if the directory does not exist, or if the
        argument refers to an HDF5 file while the code was built without HDF5 support. The HDF5
        file does not need to exist; it is created when the first bundle is written. */
    void setOutputPath(const string& value);

    /** Sets the prefix for output file names; the default is empty (i.e. no prefix). */
    void setOutputPrefix(string value);

    /** Returns the prefix for output file names. */
    string outputPrefix() const;

    /** This function determines the candidate locations for the output file with the specified
        name, as described in the class header, and registers the name of that file as being
        written on behalf of the specified simulation item. If the name has been registered
        before, the function logs a warning, as described in the class header. A class that writes
        an output file calls this function exactly once for that file, passing the simulation item
        on whose behalf the file is written; it must not change the name of the file obtained from
        the function. The function is not thread-safe; it must be called from serial code.

        The prefix returned by outputPrefix() is inserted in front of the specified name, separated
        by an underscore. The function returns a pair of strings.

        The first string is the absolute path of a regular file with that name in the output
        directory. It is always nonempty, because some output files, such as the log file, are
        written as a regular file regardless of the output path.

        The second string is the path of a bundle with that name in the HDF5 file and suite
        included in the output path, in the form <tt>\<filepath\>:\<suite\>/\<bundle\></tt>
        expected by the H5Lib::createBundle() function, or the empty string if the output path does
        not include an HDF5 file.

        The function performs no input/output operations. A caller that supports writing a bundle
        is expected to write the bundle if the second string is nonempty, and the regular file
        otherwise. */
    std::pair<string, string> output(const SimulationItem* item, const string& name);

    /** This function returns the same candidate locations for the output file with the specified
        name as the output() function, without registering the name of the file. It is intended for
        locating the log file and the parameter file, for locating a file that has already been
        written, and for determining whether the output path includes an HDF5 file, as described in
        the class header. */
    std::pair<string, string> outputLocation(const string& name) const;

    /** This function returns the directory portion of the specified input or output path, which
        can have the form <tt>\<dir\>/\<hdf\>:\<suite\></tt> described in the class header. If
        the path does not refer to an HDF5 file, the function returns the path itself. It allows
        locating the output directory before the output path has been set, for example to report an
        error that occurred while constructing the simulation. */
    static string directoryPath(const string& path);

    //========================= I/O cutover ========================

public:
    /** This function specifies whether other simulations run in parallel in the same process, as
        with the \c -s command line option of SKIRT. The default is false. The function must be
        called before the cutover. */
    void setSharedProcess(bool shared);

    /** This function records that the simulation has reached the I/O cutover, as described in the
        class header. Before doing so, if the process runs a single simulation, it verifies that the
        process holds no open HDF5 objects. If the output path includes an HDF5 file, it then
        verifies in the root process that this file can be written, creating the file if it does
        not exist. The function throws a fatal error if one of these checks fails, in which case the
        cutover is not recorded. Once the cutover has been recorded, a further call has no effect.
        */
    void cutOver();

    /** This function returns true if the simulation has reached the I/O cutover, and false
        otherwise. */
    bool hasCutOver() const;

    /** This function throws a fatal error if the simulation has reached the I/O cutover. A client
        calls it before it opens the bundle with the specified path for reading, or tests whether
        it exists. */
    void verifyBundleInput(const string& bundlePath) const;

    /** This function throws a fatal error if the simulation has not yet reached the I/O cutover. A
        client calls it before it opens the bundle with the specified path for writing. */
    void verifyBundleOutput(const string& bundlePath) const;

    //======================== Resource files =======================

public:
    /** This function returns true if a resource with the specified filename is available, and
        false otherwise. The filename should \em not include any directory segments (just the base
        filename and filename extension). The function searches the list of available resource
        files as described in the class header. */
    static bool hasResource(string name);

    /** This function returns the absolute canonical path for a resource file with the specified
        filename. The filename should \em not include any directory segments (just the base
        filename and filename extension). The function searches the list of available resource
        files as described in the class header. If the specified resource file cannot be located, a
        fatal error is thrown. */
    static string resource(const string& name);

    /** This function returns the filename (without directory segments) for a resource file with
        the specified type and with a filename including the specified segments. The function
        searches the list of available resource files as described in the class header.

        For a resource file to be considered by this function, the end of its filename (including
        the filename extension) must match the specified type string, and the filename must also
        contain each of the specified segments. If no or multiple resources files match these
        requirements, the function throws a fatal error. If a single resource file matches, the
        function returns its filename. */
    static string resourceName(const string& type, const vector<string>& segments);

    //======================== Resource packs =======================

public:
    /** This function returns a list of the names of all expected packs, in the order listed in the
        \c ExpectedResources.txt resource file, or the empty list if that resource file is missing
        or empty. */
    static vector<string> expectedPacks();

    /** This function returns the expected version number for the resource pack with the specified
        name, or zero if the name is not in the list of expected packs. */
    static int expectedPackVersion(const string& name);

    /** This function returns the version number for the installed resource pack with the specified
        name, or zero if the pack is not installed. */
    static int installedPackVersion(const string& name);

    //======================== Data Members ========================

private:
    string _inputPath;       // absolute canonical path of the input directory, with trailing slash
    string _inputHdf5Path;   // absolute path of the HDF5 file included in the input path, or empty if there is none
    string _inputSuite;      // path of the suite in that HDF5 file, without leading or trailing slash, or empty
    string _outputPath;      // absolute canonical path of the output directory, with trailing slash
    string _outputHdf5Path;  // absolute path of the HDF5 file included in the output path, or empty if there is none
    string _outputSuite;     // path of the suite in that HDF5 file, without leading or trailing slash, or empty
    string _outputPrefix;    // prefix for the names of output files
    bool _sharedProcess{false};  // true if other simulations run in parallel in the same process
    bool _cutOver{false};        // true if the simulation has reached the I/O cutover

    // an output file name registered by the output() function
    struct OutputFile
    {
        string name;                          // the name of the file, including the prefix
        const SimulationItem* item{nullptr};  // the simulation item on whose behalf the file is written
        string itemDescription;               // the type and name of that simulation item
    };
    std::unordered_map<string, OutputFile> _outputFiles;  // the registered files, indexed by lowercase name
    vector<std::pair<const SimulationItem*, const SimulationItem*>> _reportedConflicts;  // pairs of items
};

////////////////////////////////////////////////////////////////////

#endif
