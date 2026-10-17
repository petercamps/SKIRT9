/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPATIALGRIDPLOTFILE_HPP
#define SPATIALGRIDPLOTFILE_HPP

#include "Basics.hpp"
class SimulationItem;
class Units;

////////////////////////////////////////////////////////////////////

/** This class is used to write geometric information about a spatial grid to an output file in a
    format that can be easily plotted. There are two format variations for 2D and 3D information,
    respectively. The 2D format describes the intersection of a spatial grid with one of the
    coordinate planes. The 3D format fully describes all or part of the spatial cells in the grid.
    Both formats describe a sequence of polylines: the first point of a polyline is a "moveto"
    command, and each further point is a "lineto" command. All coordinates are converted to output
    units.

    The public functions of this class draw various shapes, each as one or more polylines. The
    file name passed to the constructor is passed to the FilePaths::output() function to obtain
    the candidate locations of the output file:

    - If the output path of the simulation includes an HDF5 file, the plot is written as a bundle
      in that file, as described below. The points are collected in memory, and the bundle is
      written when the file is closed.

    - Otherwise, the plot is written as a regular text file, point by point. Each line in the file
      contains two (2D) or three (3D) coordinates separated by a tab character, with 8 significant
      digits, or is empty. Consecutive nonempty lines form a polyline; an empty line ends a
      polyline, so that the next point starts a new one.

    A client must call the close() function explicitly after drawing the last shape, because
    writing a bundle may fail, and only that function can report the failure as a fatal error.
    The destructor completes a file that has not been closed, but because a destructor must not
    throw an exception, it can only log an error message, after which the simulation continues.

    In a multiprocessing environment, only the root process writes the output file; the functions
    called by other processes have no effect.

    Spatial grid plot file bundle
    -----------------------------

    A spatial grid plot file bundle holds the following datasets, where \f$N\f$ is the total
    number of points in all polylines:

    | %Dataset | Dimensions | Type | Description |
    |--|--|--|--|
    | \c points | (N, 2) or (N, 3) | 64-bit float | The coordinates of each point |
    | \c moveto | (N) | Boolean | True for the first point of a polyline, false otherwise |

    The bundle has the attribute \em unit, holding the output length unit of the coordinates, in
    addition to the standard attributes of every bundle (see the H5Lib class). A plot without any
    points has an empty \c points dataset with two columns, or with three columns if the file name
    ends with "_xyz", as for the file names used for 3D plots. All points of a plot must have the
    same number of coordinates. */
class SpatialGridPlotFile
{
    // ----------- Construction - Destruction ----------

public:
    /** The constructor determines the format of the output file, as described in the class
        header, and opens a regular output file if applicable. The arguments are (1) a simulation
        item in the hierarchy of the caller, used to retrieve the output file path, the output
        units, and an appropriate logger; and (2) the name of the file, excluding path, simulation
        prefix and filename extension. */
    SpatialGridPlotFile(const SimulationItem* item, const string& filename);

    /** The copy constructor is deleted because instances of this class cannot be copied. */
    SpatialGridPlotFile(const SpatialGridPlotFile&) = delete;

    /** The assignment operator is deleted because instances of this class cannot be copied. */
    SpatialGridPlotFile& operator=(const SpatialGridPlotFile&) = delete;

    /** This function completes the output file, if this was not already done. For a regular file,
        it closes the file. For a bundle, it writes the bundle with all points drawn so far, and
        throws a fatal error if this fails. In both cases, it logs an informational message. A
        client must call this function after drawing the last shape. */
    void close();

    /** The destructor does not throw exceptions. If the client did not call the close() function,
        the destructor calls it, unless an exception is being propagated, in which case the output
        file is left incomplete. If close() throws a fatal error, the destructor logs the error
        message instead of propagating the error. */
    ~SpatialGridPlotFile();

    // ----------- 2D functions ----------

public:
    /** This function outputs the specified 2D line segment. */
    void writeLine(double beg1, double beg2, double end1, double end2);

    /** This function outputs the 2D line segments describing the specified rectangle. */
    void writeRectangle(double min1, double min2, double max1, double max2);

    /** This function outputs a number of 2D line segments describing a circle with the specified
        radius centered at the specified coordinates. */
    void writeCircle(double x, double y, double radius);

    /** This function outputs a number of 2D line segments describing a circle centered at the
        origin with the specified radius. */
    void writeCircle(double radius);

    /** This function outputs a number of 2D line segments describing an arc centered at the
        origin with the specified radius and begin/end angles. */
    void writeArc(double radius, double begAngle, double endAngle);

    // ----------- 3D functions ----------

public:
    /** This function outputs the specified 3D line segment. */
    void writeLine(double x1, double y1, double z1, double x2, double y2, double z2);

    /** This function outputs 3D line segments describing the specified cuboid. */
    void writeCube(double x1, double y1, double z1, double x2, double y2, double z2);

    /** This function outputs a number of 3D line segments describing a circle with the specified
        radius at the specified height. */
    void writeCircle(double radius, double z);

    /** This function outputs a number of 3D line segments describing half a circle with the
        specified radius in the meridional plane with the specified angle. */
    void writeMeridionalHalfCircle(double radius, double phi);

    /** This function outputs a number of 3D line segments describing a coarse wireframe
        representation of a sphere with the specified center and radius, consisting of three
        orthogonal circles (one each in the planes parallel to the xy, xz, and yz coordinate planes
        through the center). The circles are drawn with fewer segments than writeCircle(), since
        this function is intended for plotting many small spheres rather than one accurate one. */
    void writeSphere(double x, double y, double z, double radius);

    /** This function outputs 3D line segments describing the specified polyhedron. Assuming
        the polyhedron has \f$n\f$ vertices, the first vector contains the \f$3n\f$ vertex
        coordinates as a sequence \f[x_0,y_0,z_0,x_1,y_1,z_1,...,x_{n-1},y_{n-1},y_{n-1}\f] The
        second vector contains the vertex indices for all polyhedron faces. Each vertex index is a
        number in the range \f$[0,n-1]\f$. The indices for each face are preceded by the number of
        vertices in the face. Assuming the polyhedron has \f$k\f$ faces, the index sequence looks
        like
        \f[m_0,i_{0,0},i_{0,1},...,i_{0,m_0-1},...,m_k,i_{k-1,0},i_{k-1,1},...,i_{k-1,m_k-1}\f] */
    void writePolyhedron(const vector<double>& coords, const vector<int>& indices);

    // ----------- Private helpers ----------

private:
    /** This function starts a new polyline at the specified 2D point, in output units. */
    void moveTo(double x, double y);

    /** This function continues the current polyline to the specified 2D point, in output units. */
    void lineTo(double x, double y);

    /** This function starts a new polyline at the specified 3D point, in output units. */
    void moveTo(double x, double y, double z);

    /** This function continues the current polyline to the specified 3D point, in output units. */
    void lineTo(double x, double y, double z);

    /** This function adds a point with the specified number of coordinates (2 or 3) to the output
        file, starting a new polyline if \em moveto is true. It implements the moveTo() and
        lineTo() functions. */
    void addPoint(bool moveto, int dimension, double x, double y, double z);

    /** This function returns true if this instance writes an output file, i.e. if it runs in the
        root process and the file has not yet been closed. */
    bool isActive() const { return !_isClosed && _units; }

    // ----------- Data members ----------

private:
    // a TextOutFile that offers access to its output stream, defined in the implementation file
    class PlainFile;

    // the information recorded by the constructor
    const SimulationItem* _item{nullptr};
    string _filename;
    Units* _units{nullptr};  // remains null in non-root processes
    bool _isClosed{false};

    // the regular file and the state of its current polyline, if the plot is being written as a regular file
    std::unique_ptr<PlainFile> _plain;
    bool _inPolyline{false};

    // the location of the bundle and the points collected for it, if the plot is being written as a bundle
    string _bundlePath;
    int _dimension{0};       // the number of coordinates per point, or zero if there are no points yet
    vector<double> _coords;  // the coordinates of all points, in order
    vector<char> _moveto;    // the moveto flag for each point
};

////////////////////////////////////////////////////////////////////

#endif
