# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Derive the version information from the git repository
# ------------------------------------------------------------------

# This script sets two variables describing the version of the source code, based on the release tags in the
# git repository containing this script. A release tag consists of the letter "v" followed by the version number,
# e.g. "v10.1.0"; other tags are ignored. Both lightweight and annotated tags are considered.
#
#  - BUILDINFO_PROJECT_VERSION: the most recent release tag reachable from the current commit, e.g. "v10.1.0",
#    or "unversioned" if there is no such tag or no git repository.
#  - BUILDINFO_CODE_VERSION: the output of git describe, e.g. "v10.1.0" for the tagged commit itself,
#    "v10.1.0-5-gabc1234" for commit abc1234 five commits past it, or just the abbreviated commit hash
#    "abc1234" if there is no release tag; the suffix "-dirty" indicates uncommitted changes. The value is
#    "unknown" if there is no git repository.
#
# The script can be included by a CMake list file at configure time, and by a CMake script at build time.

find_program(GIT_EXECUTABLE git)
set(BUILDINFO_PROJECT_VERSION "")
set(BUILDINFO_CODE_VERSION "")
if (GIT_EXECUTABLE)
    execute_process(COMMAND "${GIT_EXECUTABLE}" describe --tags --match "v[0-9]*" --abbrev=0
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}" OUTPUT_VARIABLE BUILDINFO_PROJECT_VERSION
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND "${GIT_EXECUTABLE}" describe --tags --match "v[0-9]*" --dirty --always
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}" OUTPUT_VARIABLE BUILDINFO_CODE_VERSION
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
if (NOT BUILDINFO_PROJECT_VERSION)
    set(BUILDINFO_PROJECT_VERSION "unversioned")
endif()
if (NOT BUILDINFO_CODE_VERSION)
    set(BUILDINFO_CODE_VERSION "unknown")
endif()

# ------------------------------------------------------------------
