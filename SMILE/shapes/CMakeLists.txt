# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds the example command-line utility 'shapes'
# ------------------------------------------------------------------

# set the target name
set(TARGET shapes)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the executable target
add_executable(${TARGET} ${SOURCES} ${HEADERS})

# add SMILE library dependencies
target_link_libraries(${TARGET} serialize schema fundamentals build)
include_directories(../serialize ../schema ../fundamentals ../build)

# adjust C++ compiler flags to our needs
include("../build/CompilerFlags.cmake")
