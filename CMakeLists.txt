# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds all targets defined in the SKIRT project
# ------------------------------------------------------------------

# ensure minimum level of functionality
cmake_minimum_required(VERSION 3.2.2...3.5 FATAL_ERROR)

# if suppported (by CMake version 3.29 and up), de-duplicate libraries on link lines based on linker capabilities
# to avoid link warnings for duplicate libraries
if(POLICY CMP0156)
  cmake_policy(SET CMP0156 NEW)
endif()

# if supported (by CMake version 3.9 and up), honor the INTERPROCEDURAL_OPTIMIZATION target property for all compilers
# rather than just for the Intel classic compiler on Linux, which is needed for link-time optimization (see below)
if(POLICY CMP0069)
  cmake_policy(SET CMP0069 NEW)
endif()

# define project
project(SKIRTproject)

# define a user-configurable option to build SKIRT
option(BUILD_SKIRT "build SKIRT, advanced radiative transfer" ON)

# define a user-configurable option to build MakeUp, which requires Qt
option(BUILD_MAKE_UP "build MakeUp, desktop GUI wizard - requires Qt5 or Qt6")

# define a user-configurable option to build doxstyle
option(BUILD_DOX_STYLE "build doxstyle, documentation block streamliner")

# define a user-configurable option to force warnings to be handled as errors
option (WARNINGS_AS_ERRORS "handle warnings as errors")

# define a user-configurable option to build with link-time optimization (LTO), which allows the compiler to optimize
# across source file and library boundaries at the cost of a longer link step; it applies to the Release build type only
option(BUILD_WITH_LTO "build with link-time optimization, if supported by the compiler" ON)

# enable link-time optimization for the Release build type of all targets, if requested and supported; this must be
# done before adding the subdirectories because a target takes on the value of this variable at the time it is created
if (BUILD_WITH_LTO AND NOT CMAKE_VERSION VERSION_LESS 3.9)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_OUTPUT LANGUAGES C CXX)
    if (APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        # GCC needs a linker plugin for LTO, which the Apple linker lacks: linking succeeds but nothing is optimized,
        # while the executable grows several-fold from carrying along the unused intermediate code
        set(LTO_SUPPORTED FALSE)
    endif()
    if (LTO_SUPPORTED)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
        message(STATUS "Link-time optimization enabled for Release builds")
    else()
        message(STATUS "Link-time optimization is not supported with this compiler and platform; building without it")
        message(VERBOSE "${LTO_OUTPUT}")
    endif()
endif()

# add all relevant subdirectories; each subdirectory defines a single target
add_subdirectory(SMILE)
if (BUILD_SKIRT)
    add_subdirectory(SKIRT)
endif()
if (BUILD_MAKE_UP)
    add_subdirectory(MakeUp)
endif()
if (BUILD_DOX_STYLE)
    add_subdirectory(Docs/doxstyle)
endif()

# set the build type to "Release" as a default; this can be overridden by the user
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "the type of build" FORCE)
endif()

# mark some obscure variables as advanced so they don't clutter the user interface
mark_as_advanced(CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT)
mark_as_advanced(CMAKE_CODEBLOCKS_EXECUTABLE CMAKE_CODEBLOCKS_MAKE_ARGUMENTS QT_QMAKE_EXECUTABLE)
mark_as_advanced(CORE_FOUNDATION)
