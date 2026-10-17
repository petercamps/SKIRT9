/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FITINSOUT_HPP
#define FITINSOUT_HPP

#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** The FITSInOut class offers static functions to read/write a 2D or 3D data stream from/to a
    standard FITS file, with support for a basic set of metadata in the header. These functions
    are private: to read a genuine input file or to write an output file in the context of a
    simulation item hierarchy, use the FITSInOutFile class, which delegates to this class. It
    locates the file relative to the input or output path of the simulation, and also supports
    files provided as a bundle in an HDF5 file. */
class FITSInOut final
{
    // allow the FITSInOutFile class to delegate to the basic read and write functions
    friend class FITSInOutFile;

public:
    // ================== Observer information ==================

    /** This basic data structure holds information on the observer associated with a data cube
        being written. These fields are included in the FITS header, or as attributes of an HDF5
        bundle, as a convenience to the user. Note: if the redshift is zero, the angular-diameter
        and luminosity distances must be equal. */
    struct ObserverInfo
    {
        double inclination{0}, azimuth{0}, roll{0};                // must be in degrees, regardless of output units
        double redshift{0};                                        // if redshift is zero, both distances must be equal
        double luminosityDistance{0}, angularDiameterDistance{0};  // distances in output units
        string distanceUnits;                                      // describes distance output units
    };

    // ================== Basic read/write ==================

private:
    /** This function reads a data frame or cube from a FITS file containing one or more data
        planes (i.e. a 2D or 3D data array). The first argument specifies a relative or absolute
        file path; a file with that name should exist. By default, the function uses the first
        nonempty data unit in the file; usually this is the primary data unit. To specify another
        data unit, add the extension name (EXTNAME) between square brackets after the file
        name. For example:

            SomePath/myFile.fits[SED]

        The subsequent function arguments serve to store data read from the file: \em data receives
        the actual values in the data cube; \em nx and \em ny receive the number of values in each
        direction, \em nz specifies the number of planes (which is equal to 1 for 2D data). The
        values in the \em data array are ordered such that the index along the x-axis varies most
        rapidly, the index along the y-axis varies less rapidly, and the index along the z-axis (if
        present) varies least rapidly. */
    static void read(const string& filepath, Array& data, int& nx, int& ny, int& nz);

    /** This function reads a single table column from a FITS file containing a two-dimensional
        table. The first argument specifies a relative or absolute file path; a file with that name
        should exist. By default, the function uses the first column of the first table in the
        file. To specify another table, add the extension name (EXTNAME) of the data unit
        containing the table between square brackets after the file name. To specify another column
        in the table, add the column name as follows:

            SomePath/myFile.fits[AXES][col time]

        The \em data argument serves to store the column data read from the file, and \em n
        specifies the number of values (i.e. rows) to be read. */
    static void readColumn(const string& filepath, Array& data, int n);

    /** This function writes a 2D data frame or a 3D data cube to the primary data unit of a FITS
        file. The x and y axes are assumed to discretize spatial extent on a regular, linear grid.
        The optional z axis often represents spectral extent and may be discretized on an arbitary
        grid. If it is present (i.e. for a 3D data cube), the function also writes the z-axis grid
        points to the FITS file as a binary table extension, so that the file is self-contained.

        The first argument specifies a relative or absolute file path; if a file with that name
        already exists, it is overwritten. The subsequent arguments specify the contents of the
        file. \em data contains the actual values in the frame or data cube; the values in this
        array must be ordered such that the index along the x-axis varies most rapidly, the index
        along the y-axis varies less rapidly, and the index along the z-axis (if present) varies
        least rapidly. \em dataUnits describes the units of the data values; \em nx and \em ny
        specify the number of values in each spatial direction, \em incx and \em incy specify the
        increment between subsequent grid points in each spatial direction, \em xc and \em yc
        specify the center of the frame(s), and \em xyUnits describes the units of the xy-grid
        increments. \em z contains the z-axis grid points (often wavelengths), and \em zUnits
        describes the units of these grid points. Finally, \em obsInfo is a pointer to an optional
        data structure holding observer information (or the null pointer). */
    static void write(const string& filepath, const Array& data, const string& dataUnits, int nx, int ny, double incx,
                      double incy, double xc, double yc, const string& xyUnits, const Array& z, const string& zUnits,
                      const ObserverInfo* obsInfo);

    /** This function writes a 2D map to the primary data unit of a FITS file. The x and y axis
        each discretizes some arbitrary quantity, for example wavelength and time lag for a
        spectral-time map. The function also writes the axis grid points to the FITS file as a
        binary table extension, so that the file is self-contained.

        The first argument specifies a relative or absolute file path; if a file with that name
        already exists, it is overwritten. The subsequent arguments specify the contents of the
        file. \em data contains the actual values in the frame; the values in this array must be
        ordered such that the index along the x-axis varies most rapidly and the index along the
        y-axis varies less rapidly. \em dataUnits describes the units of the data values. \em x and
        \em y specify the grid points in each direction, and \em xUnits and \em yUnits describe the
        units of these grid points. Finally, \em obsInfo is a pointer to an optional data structure
        holding observer information (or the null pointer). */
    static void writeMap(const string& filepath, const Array& data, const string& dataUnits, const Array& x,
                         const Array& y, const string& xUnits, const string& yUnits, const ObserverInfo* obsInfo);
};

////////////////////////////////////////////////////////////////////

#endif
