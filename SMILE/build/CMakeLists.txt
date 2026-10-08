# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Provide version and build information
# ------------------------------------------------------------------

# set the target name
set(TARGET build)

# create the library target
add_library(${TARGET} STATIC BuildInfo.cpp BuildInfo.hpp)

# add include path for Basics.hpp
include_directories(../fundamentals)

# add include path for headers created by CMake in the build directory
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# adjust C++ compiler flags to our needs
include("CompilerFlags.cmake")

# ------------------------------------------------------------------

# write a cmake script that writes the timestamp.h header containing a build timestamp and the version information
# derived from the git repository by the Version.cmake script, so that this information is current for each build
set(TS_PATH "${PROJECT_BINARY_DIR}/timestamp.cmake")
file(WRITE ${TS_PATH}  "string(TIMESTAMP BUILD_DATE \"%d/%m/%Y\")\n")
file(APPEND ${TS_PATH} "string(TIMESTAMP BUILD_TIME \"%H:%M:%S\")\n")
file(APPEND ${TS_PATH} "include(\"${CMAKE_CURRENT_SOURCE_DIR}/Version.cmake\")\n")
file(APPEND ${TS_PATH} "file(WRITE timestamp.h \"#ifndef TIMESTAMP_H\\n\")\n")
file(APPEND ${TS_PATH} "file(APPEND timestamp.h \"#define TIMESTAMP_H\\n\")\n")
file(APPEND ${TS_PATH} "file(APPEND timestamp.h \"#define BUILD_DATE \\\"\${BUILD_DATE}\\\"\\n\")\n")
file(APPEND ${TS_PATH} "file(APPEND timestamp.h \"#define BUILD_TIME \\\"\${BUILD_TIME}\\\"\\n\")\n")
file(APPEND ${TS_PATH} "file(APPEND timestamp.h \"#define PROJECT_VERSION \\\"\${BUILDINFO_PROJECT_VERSION}\\\"\\n\")\n")
file(APPEND ${TS_PATH} "file(APPEND timestamp.h \"#define CODE_VERSION \\\"\${BUILDINFO_CODE_VERSION}\\\"\\n\")\n")
file(APPEND ${TS_PATH} "file(APPEND timestamp.h \"#endif\\n\")\n")

# provide a custom target that will execute the above cmake script at build time
add_custom_target(timestamp
    COMMAND ${CMAKE_COMMAND} -P ${PROJECT_BINARY_DIR}/timestamp.cmake
    ADD_DEPENDENCIES ${PROJECT_BINARY_DIR}/timestamp.cmake)

# add a dependency to ensure that the above custom target is executed for each build
add_dependencies(${TARGET} timestamp)

# ------------------------------------------------------------------
