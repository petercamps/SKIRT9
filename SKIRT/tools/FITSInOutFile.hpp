/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FITSINOUTFILE_HPP
#define FITSINOUTFILE_HPP

#include "FITSInOut.hpp"
class SimulationItem;

////////////////////////////////////////////////////////////////////

/** The FITSInOutFile class offers static functions to read a 2D or 3D data array from an input
    file, and to write a 2D data frame or 3D data cube to an output file, where each file is
    provided either as a standard FITS file or as a bundle in an HDF5 file. The functions are
    intended for genuine input and output files in the context of a simulation item hierarchy.
    They wrap the corresponding functions of the FITSInOut class, which accesses FITS files only.

    <b>Locating the input file</b>

    The input file is located through the FilePaths::input() function, which combines the
    specified file name with the input path of the simulation. If that function returns a regular
    file candidate and the file exists, the data is read from that FITS file. Otherwise, the data
    is read from the HDF5 bundle candidate, if there is one. This candidate exists if the input
    path of the simulation includes an HDF5 file, or if the specified file name embeds a reference
    to a bundle in an HDF5 file. See the FilePaths class for more information. If the input file
    is found at neither of two candidate locations, the fatal error that reports this mentions
    both the regular file and the bundle.

    <b>Locating the output file</b>

    The output file is located through the FilePaths::output() function, which combines the
    specified file name with the output path and prefix of the simulation. If the output path of
    the simulation includes an HDF5 file, the frame or data cube is written as a bundle in that
    file, as described below. Otherwise, it is written as a regular FITS file in the output
    directory, exactly as the FITSInOut class writes it. In a multiprocessing environment, only the
    root process writes the output file.

    <b>FITS file bundle for input</b>

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

    Any other datasets in the bundle and all attributes are ignored. The bundle written for a FITS
    output file, described below, holds its frame or data cube in a dataset with the same name and
    layout. Such an output bundle can therefore be read as input, just as a FITS file written by a
    previous simulation can.

    The name of a FITS file can be followed by the name of an extension between square brackets
    to select a data unit other than the first one (see the FITSInOut class). There is no
    equivalent for a bundle: each data unit should be placed in its own bundle.

    <b>FITS file bundle for output</b>

    The bundle written for a FITS output file holds the frame or data cube together with the grid
    points of its axes, in the following datasets:

    | %Dataset | Dimensions | Type | Description |
    |--|--|--|--|
    | \c data | (ny, nx) or (nz, ny, nx) | 64-bit float | Pixel or voxel values |
    | \c x | (nx) | 64-bit float | Grid points of the x-axis |
    | \c y | (ny) | 64-bit float | Grid points of the y-axis |
    | \c z | (nz) | 64-bit float | Grid points of the z-axis, if there is one |

    Each dataset has the attribute \em unit, holding the unit string specified for the data or for
    the axis, as recorded in the corresponding card of a FITS file (BUNIT, CUNIT1, CUNIT2, and
    CUNIT3). The data values are stored as 64-bit floating point numbers, while a FITS file holds
    32-bit floating point numbers. The grid points of a regularly spaced x-axis or y-axis, which a
    FITS file describes through its CRPIX, CRVAL, and CDELT cards, are the coordinates of the
    pixel centers, \f$x_i = x_\mathrm{c} + (i + 1 - (n_x + 1)/2)\,\Delta x\f$ for \f$i = 0, \dots,
    n_x-1\f$. The grid points of an irregular axis, which a FITS file lists in a table extension,
    are stored as specified. A frame without z-axis has no \c z dataset; a data cube with a single
    frame, which a FITS file holds as a 3D array, has two-dimensional \c data and a \c z dataset
    with a single grid point.

    The bundle has the attribute \em description, holding the description of the data that is
    also used in the log message, in addition to the standard attributes of every bundle (see the
    H5Lib class). If observer information is provided, the bundle also has the 64-bit floating
    point attributes \em inclination, \em azimuth, and \em roll (in degrees), \em redshift,
    \em luminosity_distance, and \em angular_diameter_distance, and the string attribute
    \em distance_unit, with the values of the corresponding fields of the ObserverInfo structure.
    */
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

    /** This function writes a 2D data frame or a 3D data cube to an output file in the context of
        the simulation item hierarchy specified through the first argument. This allows the
        function to issue log messages using the human-readable description of the data given as
        the second argument, to skip writing in non-root processes of a multi-process simulation,
        and to locate the output file from the file name specified as the third argument, as
        described in the class header. The output filename should \em not include the filename
        extension nor the simulation prefix.

        The remaining arguments specify the contents of the file. \em data contains the actual
        values in the frame or data cube; the values in this array must be ordered such that the
        index along the x-axis varies most rapidly, the index along the y-axis varies less
        rapidly, and the index along the z-axis (if present) varies least rapidly. \em dataUnits
        describes the units of the data values; \em nx and \em ny specify the number of values in
        each spatial direction, \em incx and \em incy specify the increment between subsequent
        grid points in each spatial direction, \em xc and \em yc specify the center of the
        frame(s), and \em xyUnits describes the units of the xy-grid increments. \em z contains
        the z-axis grid points (often wavelengths), and \em zUnits describes the units of these
        grid points. If \em z is empty, the data is a single frame; otherwise it is a data cube
        with a frame for each grid point. Finally, \em obsInfo is a pointer to an optional data
        structure holding observer information (or the null pointer). The arguments describing
        the z-axis and observer info may be omitted when writing a 2D data frame. */
    static void write(const SimulationItem* item, const string& description, const string& filename, const Array& data,
                      const string& dataUnits, int nx, int ny, double incx, double incy, double xc, double yc,
                      const string& xyUnits, const Array& z = Array(), const string& zUnits = string(),
                      const FITSInOut::ObserverInfo* obsInfo = nullptr);

    /** This function writes a 2D map to an output file in the context of the simulation item
        hierarchy specified through the first argument. The first three arguments have the same
        meaning as for the write() function. The x and y axes each discretize some arbitrary
        quantity, for example wavelength and time lag for a spectral-time map. \em data contains
        the actual values in the map; the values in this array must be ordered such that the index
        along the x-axis varies most rapidly and the index along the y-axis varies less rapidly.
        \em dataUnits describes the units of the data values. \em x and \em y specify the grid
        points in each direction, and \em xUnits and \em yUnits describe the units of these grid
        points. Finally, \em obsInfo is a pointer to an optional data structure holding observer
        information (or the null pointer). */
    static void writeMap(const SimulationItem* item, const string& description, const string& filename,
                         const Array& data, const string& dataUnits, const Array& x, const Array& y,
                         const string& xUnits, const string& yUnits, const FITSInOut::ObserverInfo* obsInfo = nullptr);
};

////////////////////////////////////////////////////////////////////

#endif
