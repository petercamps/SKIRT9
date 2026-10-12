# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library with a thin C++ wrapper around the HDF5 C library
# if requested and available, or a stub that reports a fatal error
# when attempting to open an HDF5 file if not. Concentrating all
# HDF5 calls in this library allows all other code to be compiled
# without any dependency on HDF5.
# ------------------------------------------------------------------

# set the target name; it must differ from "hdf5", which is the name of a target exported by the HDF5 library itself
set(TARGET h5)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} fundamentals)
include_directories(../../SMILE/fundamentals)

# define a user-configurable option to build with HDF5 support,
# which requires the HDF5 C library to be installed on the system
option(BUILD_HDF5 "build with HDF5 support - requires the HDF5 C library")

# if requested, find the HDF5 C library and configure accordingly
if (BUILD_HDF5)
    # find and configure HDF5; the headers are included as system headers to avoid warnings caused by their contents
    find_package(HDF5 REQUIRED COMPONENTS C)
    target_link_libraries(${TARGET} ${HDF5_C_LIBRARIES})
    include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
    # add the C++ preprocessor symbol that triggers compilation with HDF5, and any symbols required by the library
    add_definitions(-DBUILD_HDF5 ${HDF5_DEFINITIONS})
endif()

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
