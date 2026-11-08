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
    /** Returns the project version, i.e. the most recent release tag reachable from the commit
        from which this executable was built (e.g. "v10.1.0"), or "unversioned" if there is no
        such tag or no git repository. */
    static string projectVersion();

    /** Returns the full version of the code from which this executable was built, i.e. the
        output of git describe. This is the release tag for a build of the tagged commit itself
        (e.g. "v10.1.0"), or the release tag followed by the number of commits past it and the
        abbreviated commit hash (e.g. "v10.1.0-5-gabc1234"). If there is no release tag, the
        abbreviated commit hash follows "unversioned-g" (e.g. "unversioned-gabc1234"). The suffix
        "-dirty" indicates uncommitted changes, and "unversioned" means that there is no git
        repository. In all cases, the full version starts with the project version. */
    static string codeVersion();

    /** Returns the time of the current build as a string formatted for human consumption (e.g.
        "built on 07/11/2026 at 21:50:38"). */
    static string timestamp();

    /** Returns a description of the current build for human consumption, consisting of the full
        version followed by the build time in parentheses (e.g. "v10.1.0-5-gabc1234 (built on
        07/11/2026 at 21:50:38)"). Callers usually precede it by the name of the program. */
    static string description();
};

////////////////////////////////////////////////////////////////////

#endif
