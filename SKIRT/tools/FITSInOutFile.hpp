/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FITSINOUTFILE_HPP
#define FITSINOUTFILE_HPP

#include "Array.hpp"
class SimulationItem;

////////////////////////////////////////////////////////////////////

/** The FITSInOutFile class offers a static function to read a 2D or 3D data array from an input
    file that is provided either as a standard FITS file or as a bundle in an HDF5 file. It is
    intended for reading genuine input files in the context of a simulation item hierarchy. The
    function wraps the corresponding function of the FITSInOut class, which accesses FITS files
    only.

    <b>Locating the input file</b>

    The input file is located through the FilePaths::input() function, which combines the
    specified file name with the input path of the simulation. If that function returns a regular
    file candidate and the file exists, the data is read from that FITS file. Otherwise, the data
    is read from the HDF5 bundle candidate, if there is one. This candidate exists if the input
    path of the simulation includes an HDF5 file, or if the specified file name embeds a reference
    to a bundle in an HDF5 file. See the FilePaths class for more information. If the input file
    is found at neither of two candidate locations, the fatal error that reports this mentions
    both the regular file and the bundle.

    <b>FITS file bundle</b>

    The bundle that corresponds to a FITS input file is an HDF5 group that contains the following
    dataset:

    | %Dataset | Dimensions | Type | Description |
    |--|--|--|--|
    | \c data | (ny, nx) or (nz, ny, nx) | 64-bit float | Pixel or voxel values |

    Just as for a FITS file, the values are stored such that the index along the x-axis varies
    most rapidly, the index along the y-axis varies less rapidly, and the index along the z-axis
    (if present) varies least rapidly. The dataset can be stored in any integer or floating point
    type that can be converted to a 64-bit floating point number without loss of information, and
    it can be chunked or compressed; the HDF5 library handles this transparently. The bundle needs
    no attributes. None of the metadata that a FITS file might carry (such as the pixel scale) is
    read from either type of file; this information is configured in the ski file.

    Any other datasets in the bundle and all attributes are ignored. The bundle that corresponds
    to a FITS output file holds its frame or data cube in a dataset with the same name and layout,
    next to datasets with the coordinates of its axes. Such an output bundle can therefore be read
    as input, just as a FITS file written by a previous simulation can.

    The name of a FITS file can be followed by the name of an extension between square brackets
    to select a data unit other than the first one (see the FITSInOut class). There is no
    equivalent for a bundle: each data unit should be placed in its own bundle. */
class FITSInOutFile final
{
public:
    /** This function reads a 2D or 3D data array from an input file in the context of the
        simulation item hierarchy specified through the first argument. This allows the function
        to issue log messages and to locate the input file from the file name specified as the
        second argument, as described in the class header. The file name should include the
        filename extension, if any.

        The remaining arguments receive the data read from the file: \em data receives the values
        in the data frame or cube; \em nx and \em ny receive the number of values in each
        direction, and \em nz receives the number of planes (which is equal to 1 for 2D data). The
        values in the \em data array are ordered such that the index along the x-axis varies most
        rapidly, the index along the y-axis varies less rapidly, and the index along the z-axis
        (if present) varies least rapidly. */
    static void read(const SimulationItem* item, const string& filename, Array& data, int& nx, int& ny, int& nz);
};

////////////////////////////////////////////////////////////////////

#endif
