/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BUILDINFO_HPP
#define BUILDINFO_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** This class provides information about the current build, such as the build time and the
    version of the source code. This information is provided to the source code at build time by
    CMake. The version information is derived from the release tags in the git repository, which
    consist of the letter "v" followed by the version number (e.g. "v10.1.0"). */
class BuildInfo final
{
public:
    /** Returns the time of current build as a string formatted for human consumption. */
    static string timestamp();

    /** Returns the project version, i.e. the most recent release tag reachable from the commit
        from which this executable was built (e.g. "v10.1.0"), or "unversioned" if there is no
        such tag or no git repository. */
    static string projectVersion();

    /** Returns a description of the code version from which this executable was built, i.e.
        "git" followed by the output of git describe. This is the release tag for a build of the
        tagged commit itself (e.g. "git v10.1.0"), the release tag followed by the number of
        commits past it and the abbreviated commit hash (e.g. "git v10.1.0-5-gabc1234"), or the
        abbreviated commit hash if there is no release tag (e.g. "git abc1234"). The suffix
        "-dirty" indicates uncommitted changes, and "git unknown" means that there is no git
        repository. */
    static string codeVersion();
};

////////////////////////////////////////////////////////////////////

#endif
