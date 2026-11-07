/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "HEALPixSkyInstrument.hpp"
#include "FatalError.hpp"
#include "FluxRecorder.hpp"
#include "Log.hpp"

////////////////////////////////////////////////////////////////////

void HEALPixSkyInstrument::setupSelfBefore()
{
    PeelOffInstrument::setupSelfBefore();

    // verify consistency of vector properties
    Vec co(_Ox - _Cx, _Oy - _Cy, _Oz - _Cz);  // vector in direction from crosshair to observer
    Vec up(_Ux, _Uy, _Uz);                    // vector in the upward direction
    if (co.norm() < 1e-20) throw FATALERROR("Crosshair is too close to observer");
    if (up.norm() < 1e-20) throw FATALERROR("Upwards direction cannot be null vector");
    if (Vec::cross(co, up).norm() < 1e-20)
        throw FATALERROR("Upwards direction cannot be parallel to viewing direction");

    // create the HEALPix tessellation
    _tessellation = HEALPixTessellation(_order);

    // we map the 12 * Nside^2 HEALPix pixels to an almost square image with 4 * Nside * (4 * Nside - 1) pixels
    _Nx = _tessellation.maxPixelsInRing();
    _Ny = _tessellation.numRings();

    // setup the transformation from world to observer coordinates

    // translate to observer position
    _transform.translate(-_Ox, -_Oy, -_Oz);

    // unit vector in direction from crosshair to observer
    Vec kn(_Ox - _Cx, _Oy - _Cy, _Oz - _Cz);
    kn /= kn.norm();
    double a = kn.x();
    double b = kn.y();
    double c = kn.z();

    // rotate from world to observer coordinates just as for the perspective transformation
    double v = sqrt(b * b + c * c);
    if (v > 0.3)
    {
        _transform.rotateX(c / v, -b / v);
        _transform.rotateY(v, -a);
        double k = (b * b + c * c) * _Ux - a * b * _Uy - a * c * _Uz;
        double l = c * _Uy - b * _Uz;
        double u = sqrt(k * k + l * l);
        _transform.rotateZ(l / u, -k / u);
    }
    else
    {
        v = sqrt(a * a + c * c);
        _transform.rotateY(c / v, -a / v);
        _transform.rotateX(v, -b);
        double k = c * _Ux - a * _Uz;
        double l = (a * a + c * c) * _Uy - a * b * _Ux - b * c * _Uz;
        double u = sqrt(k * k + l * l);
        _transform.rotateZ(l / u, -k / u);
    }

    // rotate the axes into the alignment appropriate for our purposes (z-axis up, x-axis towards crosshair)
    _transform.rotateX(0., -1.);
    _transform.rotateZ(0., 1.);

    // determine the solid angle corresponding to each pixel
    // each HEALPix pixel has the same spherical area:  4 pi/(12 Nside^2)
    double omega = _tessellation.pixelSolidAngle();

    // determine the (approximate) scale and center of the output map axes
    // x: theta in [0, pi]
    // y: phi in [0, 2*pi]
    double incx = 2. * M_PI / _Nx;
    constexpr double ctrx = M_PI;
    double incy = M_PI / _Ny;
    constexpr double ctry = M_PI / 2.;

    // configure flux recorder
    instrumentFluxRecorder()->includeSurfaceBrightnessForLocal(_Nx, _Ny, omega, incx, incy, ctrx, ctry, "posangle");
}

////////////////////////////////////////////////////////////////////

bool HEALPixSkyInstrument::hasSameSightLine(const PeelOffInstrument* other) const
{
    auto sky = dynamic_cast<const HEALPixSkyInstrument*>(other);
    return sky && radius() == sky->radius() && observerX() == sky->observerX() && observerY() == sky->observerY()
           && observerZ() == sky->observerZ() && crossX() == sky->crossX() && crossY() == sky->crossY()
           && crossZ() == sky->crossZ() && upX() == sky->upX() && upY() == sky->upY() && upZ() == sky->upZ();
}

////////////////////////////////////////////////////////////////////

Direction HEALPixSkyInstrument::bfkobs(Position bfr) const
{
    // vector and distance from launch to observer
    Vec k = Vec(_Ox, _Oy, _Oz) - bfr;
    double d = k.norm();

    // if the distance is very small, return something arbitrary - the photon packet will not be detected anyway
    if (d < _radius) return Direction();

    // otherwise return a unit vector in the direction from launch to observer
    return Direction(k / d, false);
}

////////////////////////////////////////////////////////////////////

Direction HEALPixSkyInstrument::bfky(Position bfr) const
{
    // vector and distance from launch to observer
    Vec k = Vec(_Ox, _Oy, _Oz) - bfr;
    double d = k.norm();

    // if the distance is very small, return something arbitrary - the photon packet will not be detected anyway
    if (d < _radius) return Direction();

    // vector in the plane normal to the line launch-observer
    // oriented along the projection of the up direction in that plane
    Vec ku(_Ux, _Uy, _Uz);

    // return unit vector along y-axis
    return Direction(Vec::cross(k, Vec::cross(ku, k)), true);
}

////////////////////////////////////////////////////////////////////

PeelOffInstrument::Detection HEALPixSkyInstrument::locate(Position bfr) const
{
    // transform launch position from world to observer coordinates
    Position p(_transform.transform(bfr));

    // get the spherical coordinates of the launch position relative to the observer
    double d, theta, phi;
    p.spherical(d, theta, phi);

    // if the radial distance is very small, ignore the photon packet
    Detection detection;
    if (d < _radius) return detection;

    // p.spherical() returns theta in [0, pi] and phi in [-pi, pi]
    // the HEALPix mapping algorithm expects theta in [0, pi] and phi in [0, 2*pi]
    // we center the output image on the crosshair direction by offsetting phi
    phi += M_PI;

    // determine the ring index (vertical pixel index) and the pixel-in-ring index (horizontal pixel index)
    auto [i, j] = _tessellation.locate(theta, phi);

    // detect the photon packet in the appropriate pixel of the data cube and at the appropriate distance
    detection.pixel = i + _Nx * j;
    detection.distance = d;
    return detection;
}

////////////////////////////////////////////////////////////////////
