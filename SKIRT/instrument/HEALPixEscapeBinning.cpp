/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "HEALPixEscapeBinning.hpp"
#include "ColumnOutFile.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void HEALPixEscapeBinning::setupSelfBefore()
{
    EscapeBinning::setupSelfBefore();

    _tessellation = HEALPixTessellation(_order);
}

////////////////////////////////////////////////////////////////////

int HEALPixEscapeBinning::numBins() const
{
    return _tessellation.numPixels();
}

////////////////////////////////////////////////////////////////////

int HEALPixEscapeBinning::bin(Direction bfk) const
{
    // get the spherical coordinates of the direction, with the azimuth in [0, 2 pi]
    double theta = acos(std::clamp(bfk.kz(), -1., 1.));
    double phi = atan2(bfk.ky(), bfk.kx());
    if (phi < 0.) phi += 2. * M_PI;

    // locate the pixel and return its index in the ring ordering scheme
    auto [i, j] = _tessellation.locate(theta, phi);
    return _tessellation.ringPixelIndex(i, j);
}

////////////////////////////////////////////////////////////////////

double HEALPixEscapeBinning::skyFraction(int /*bin*/) const
{
    return 1. / _tessellation.numPixels();
}

////////////////////////////////////////////////////////////////////

string HEALPixEscapeBinning::describeBin(int bin) const
{
    auto units = find<Units>();
    auto [i, j] = _tessellation.ringPixelLocation(bin);
    auto [theta, phi] = _tessellation.pixelCenter(i, j);
    return "HEALPix pixel in ring " + std::to_string(j) + " at position " + std::to_string(i)
           + ", centered at inclination " + StringUtils::toString(units->oposangle(theta), 'g', 4) + " "
           + units->uposangle() + ", azimuth " + StringUtils::toString(units->oposangle(phi), 'g', 4) + " "
           + units->uposangle();
}

////////////////////////////////////////////////////////////////////

void HEALPixEscapeBinning::addBinTableColumns(ColumnOutFile& file) const
{
    auto units = find<Units>();
    file.addColumn("ring index", "", 'd');
    file.addColumn("pixel-in-ring index", "", 'd');
    file.addColumn("center inclination", units->uposangle());
    file.addColumn("center azimuth", units->uposangle());
}

////////////////////////////////////////////////////////////////////

vector<double> HEALPixEscapeBinning::binTableValues(int bin) const
{
    auto units = find<Units>();
    auto [i, j] = _tessellation.ringPixelLocation(bin);
    auto [theta, phi] = _tessellation.pixelCenter(i, j);
    return {static_cast<double>(j), static_cast<double>(i), units->oposangle(theta), units->oposangle(phi)};
}

////////////////////////////////////////////////////////////////////

EscapeBinning::FrameLayout HEALPixEscapeBinning::frameLayout() const
{
    // use the same layout and axes as the HEALPixSkyInstrument, without its offset for the crosshair direction:
    // x: phi in [0, 2*pi], y: theta in [0, pi]
    FrameLayout layout;
    layout.numPixelsX = _tessellation.maxPixelsInRing();
    layout.numPixelsY = _tessellation.numRings();
    layout.incrementX = 2. * M_PI / layout.numPixelsX;
    layout.incrementY = M_PI / layout.numPixelsY;
    layout.centerX = M_PI;
    layout.centerY = M_PI / 2.;
    layout.quantityXY = "posangle";
    return layout;
}

////////////////////////////////////////////////////////////////////

int HEALPixEscapeBinning::framePixel(int bin) const
{
    auto [i, j] = _tessellation.ringPixelLocation(bin);
    return i + _tessellation.maxPixelsInRing() * j;
}

////////////////////////////////////////////////////////////////////
