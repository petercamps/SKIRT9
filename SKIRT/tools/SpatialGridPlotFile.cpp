/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpatialGridPlotFile.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "H5Lib.hpp"
#include "Log.hpp"
#include "ProcessManager.hpp"
#include "SimulationItem.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"
#include <exception>
#include <iomanip>

////////////////////////////////////////////////////////////////////

// a TextOutFile that offers access to its output stream, so that points can be written without formatting overhead
class SpatialGridPlotFile::PlainFile : public TextOutFile
{
public:
    using TextOutFile::TextOutFile;
    std::ofstream& out() { return _out; }
};

////////////////////////////////////////////////////////////////////

SpatialGridPlotFile::SpatialGridPlotFile(const SimulationItem* item, const string& filename)
    : _item(item), _filename(filename)
{
    // only the root process writes the output file
    if (!ProcessManager::isRoot()) return;
    _units = item->find<Units>();

    // determine the candidate locations for the output file: a regular file and a bundle in an HDF5 file
    auto bundlePath = item->find<FilePaths>()->output(filename + ".dat").second;

    // open a regular file, or remember the location of the bundle
    if (bundlePath.empty())
    {
        _plain = std::make_unique<PlainFile>(item, filename, "data to plot the spatial grid");
        _plain->out() << std::setprecision(8);
    }
    else
    {
        _bundlePath = bundlePath;
    }
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::close()
{
    if (!isActive()) return;
    _isClosed = true;

    // end the last polyline and close the regular file, which logs its own message
    if (_plain)
    {
        if (_inPolyline) _plain->out() << '\n';
        _plain->close();
        _plain.reset();
    }

    // or write the bundle with the collected points
    else
    {
        size_t dimension =
            static_cast<size_t>(_dimension ? _dimension : (StringUtils::endsWith(_filename, "_xyz") ? 3 : 2));
        size_t numPoints = _moveto.size();

        H5BundleW bundle = H5Lib::createBundle(_bundlePath);
        bundle.setStringAttribute("unit", _units->ulength());
        Array points(_coords.data(), _coords.size());
        vector<double>().swap(_coords);
        bundle.createDoubleDataset("points", points, {numPoints, dimension});
        bundle.createBoolDataset("moveto", _moveto);
        vector<char>().swap(_moveto);

        _item->find<Log>()->info(_item->typeAndName() + " wrote data to plot the spatial grid to bundle "
                                 + _bundlePath);
    }
}

////////////////////////////////////////////////////////////////////

SpatialGridPlotFile::~SpatialGridPlotFile()
{
    // complete the output file if the client did not call close(), unless an exception is being propagated
    if (isActive() && std::uncaught_exceptions() == 0)
    {
        // a destructor must not throw, so log the message of a fatal error and ignore other exceptions
        try
        {
            try
            {
                close();
            }
            catch (FatalError& error)
            {
                auto log = _item->find<Log>();
                for (const string& line : error.message()) log->error(line);
            }
        }
        catch (...)
        {}
    }
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeLine(double beg1, double beg2, double end1, double end2)
{
    if (!isActive()) return;

    beg1 = _units->olength(beg1);
    beg2 = _units->olength(beg2);
    end1 = _units->olength(end1);
    end2 = _units->olength(end2);
    moveTo(beg1, beg2);
    lineTo(end1, end2);
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeRectangle(double min1, double min2, double max1, double max2)
{
    if (!isActive()) return;

    min1 = _units->olength(min1);
    min2 = _units->olength(min2);
    max1 = _units->olength(max1);
    max2 = _units->olength(max2);
    moveTo(min1, min2);
    lineTo(min1, max2);
    lineTo(max1, max2);
    lineTo(max1, min2);
    lineTo(min1, min2);
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeCircle(double x, double y, double radius)
{
    if (!isActive()) return;

    x = _units->olength(x);
    y = _units->olength(y);
    radius = _units->olength(radius);

    for (int l = 0; l <= 360; l++)
    {
        double phi = l * M_PI / 180;
        addPoint(l == 0, 2, x + radius * cos(phi), y + radius * sin(phi), 0.);
    }
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeCircle(double radius)
{
    if (!isActive()) return;

    radius = _units->olength(radius);

    for (int l = 0; l <= 360; l++)
    {
        double phi = l * M_PI / 180;
        addPoint(l == 0, 2, radius * cos(phi), radius * sin(phi), 0.);
    }
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeArc(double radius, double begAngle, double endAngle)
{
    if (!isActive()) return;

    radius = _units->olength(radius);

    constexpr double dphi = M_PI / 180.;
    bool first = true;
    for (double phi = begAngle; phi < endAngle; phi += dphi)
    {
        addPoint(first, 2, radius * cos(phi), radius * sin(phi), 0.);
        first = false;
    }
    addPoint(first, 2, radius * cos(endAngle), radius * sin(endAngle), 0.);
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeLine(double x1, double y1, double z1, double x2, double y2, double z2)
{
    if (!isActive()) return;

    x1 = _units->olength(x1);
    y1 = _units->olength(y1);
    z1 = _units->olength(z1);
    x2 = _units->olength(x2);
    y2 = _units->olength(y2);
    z2 = _units->olength(z2);
    moveTo(x1, y1, z1);
    lineTo(x2, y2, z2);
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeCube(double x1, double y1, double z1, double x2, double y2, double z2)
{
    if (!isActive()) return;

    x1 = _units->olength(x1);
    y1 = _units->olength(y1);
    z1 = _units->olength(z1);
    x2 = _units->olength(x2);
    y2 = _units->olength(y2);
    z2 = _units->olength(z2);
    moveTo(x1, y1, z1);
    lineTo(x2, y1, z1);
    lineTo(x2, y2, z1);
    lineTo(x1, y2, z1);
    lineTo(x1, y1, z1);

    moveTo(x1, y1, z2);
    lineTo(x2, y1, z2);
    lineTo(x2, y2, z2);
    lineTo(x1, y2, z2);
    lineTo(x1, y1, z2);

    moveTo(x1, y1, z1);
    lineTo(x1, y1, z2);
    moveTo(x2, y1, z1);
    lineTo(x2, y1, z2);
    moveTo(x2, y2, z1);
    lineTo(x2, y2, z2);
    moveTo(x1, y2, z1);
    lineTo(x1, y2, z2);
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeCircle(double radius, double z)
{
    if (!isActive()) return;

    radius = _units->olength(radius);
    z = _units->olength(z);

    for (int l = 0; l <= 360; l++)
    {
        double phi = l * M_PI / 180;
        addPoint(l == 0, 3, radius * cos(phi), radius * sin(phi), z);
    }
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeMeridionalHalfCircle(double radius, double phi)
{
    if (!isActive()) return;

    radius = _units->olength(radius);

    for (int l = 0; l <= 180; l++)
    {
        double theta = l * M_PI / 180;
        double x = radius * cos(phi) * sin(theta);
        double y = radius * sin(phi) * sin(theta);
        double z = radius * cos(theta);
        addPoint(l == 0, 3, x, y, z);
    }
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writeSphere(double x, double y, double z, double radius)
{
    if (!isActive()) return;

    x = _units->olength(x);
    y = _units->olength(y);
    z = _units->olength(z);
    radius = _units->olength(radius);

    constexpr int N = 24;  // segments per circle; coarse, since there may be many spheres to plot

    // circle in the plane parallel to the xy plane
    for (int l = 0; l <= N; l++)
    {
        double phi = l * 2. * M_PI / N;
        addPoint(l == 0, 3, x + radius * cos(phi), y + radius * sin(phi), z);
    }

    // circle in the plane parallel to the xz plane
    for (int l = 0; l <= N; l++)
    {
        double phi = l * 2. * M_PI / N;
        addPoint(l == 0, 3, x + radius * cos(phi), y, z + radius * sin(phi));
    }

    // circle in the plane parallel to the yz plane
    for (int l = 0; l <= N; l++)
    {
        double phi = l * 2. * M_PI / N;
        addPoint(l == 0, 3, x, y + radius * cos(phi), z + radius * sin(phi));
    }
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::writePolyhedron(const vector<double>& coords, const vector<int>& indices)
{
    if (!isActive()) return;

    unsigned int k = 0;
    while (k < indices.size())
    {
        int numvertices = indices[k++];
        int firstindex = indices[k];
        for (int i = 0; i < numvertices; i++)
        {
            int currentindex = indices[k++];
            double x = _units->olength(coords[3 * currentindex]);
            double y = _units->olength(coords[3 * currentindex + 1]);
            double z = _units->olength(coords[3 * currentindex + 2]);
            addPoint(i == 0, 3, x, y, z);
        }
        double x = _units->olength(coords[3 * firstindex]);
        double y = _units->olength(coords[3 * firstindex + 1]);
        double z = _units->olength(coords[3 * firstindex + 2]);
        lineTo(x, y, z);
    }
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::moveTo(double x, double y)
{
    addPoint(true, 2, x, y, 0.);
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::lineTo(double x, double y)
{
    addPoint(false, 2, x, y, 0.);
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::moveTo(double x, double y, double z)
{
    addPoint(true, 3, x, y, z);
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::lineTo(double x, double y, double z)
{
    addPoint(false, 3, x, y, z);
}

////////////////////////////////////////////////////////////////////

void SpatialGridPlotFile::addPoint(bool moveto, int dimension, double x, double y, double z)
{
    // for a regular file, end the current polyline with an empty line before starting a new one
    if (_plain)
    {
        auto& out = _plain->out();
        if (moveto && _inPolyline) out << '\n';
        out << x << '\t' << y;
        if (dimension == 3) out << '\t' << z;
        out << '\n';
        _inPolyline = true;
    }

    // for a bundle, collect the point; the first point always starts a polyline
    else
    {
        if (!_dimension) _dimension = dimension;
        if (dimension != _dimension)
            throw FATALERROR("Points with two and three coordinates are mixed in spatial grid plot " + _filename);
        _coords.push_back(x);
        _coords.push_back(y);
        if (dimension == 3) _coords.push_back(z);
        _moveto.push_back(moveto || _moveto.empty());
    }
}

////////////////////////////////////////////////////////////////////
