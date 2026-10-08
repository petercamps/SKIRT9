/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PerspectiveInstrument.hpp"
#include "FatalError.hpp"
#include "FluxRecorder.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    double norm(double x, double y, double z)
    {
        return sqrt(x * x + y * y + z * z);
    }
}

////////////////////////////////////////////////////////////////////

void PerspectiveInstrument::setupSelfBefore()
{
    Instrument::setupSelfBefore();

    // verify attribute values
    if (_Ux == 0 && _Uy == 0 && _Uz == 0) throw FATALERROR("Upwards direction cannot be null vector");

    // unit vector in direction from crosshair position to viewport origin
    double G = norm(_Vx - _Cx, _Vy - _Cy, _Vz - _Cz);
    if (G < 1e-20) throw FATALERROR("Crosshair is too close to viewport origin");
    double a = (_Vx - _Cx) / G;
    double b = (_Vy - _Cy) / G;
    double c = (_Vz - _Cz) / G;

    // pixel width and height (assuming square pixels)
    _s = _Sx / _Nx;

    // eye position
    _Ex = _Vx + _Fe * a;
    _Ey = _Vy + _Fe * b;
    _Ez = _Vz + _Fe * c;

    // unit vectors along viewport's x and y axes
    Vec kn(_Vx - _Cx, _Vy - _Cy, _Vz - _Cz);
    Vec ku(_Ux, _Uy, _Uz);
    Vec ky = Vec::cross(kn, Vec::cross(ku, kn));
    _bfky = Direction(ky, true);

    // the perspective transformation

    // from world to eye coordinates
    _transform.translate(-_Ex, -_Ey, -_Ez);
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
    _transform.scale(1., 1., -1.);

    // from eye to viewport coordinates
    _transform.perspectiveZ(_Fe);

    // from viewport to pixel coordinates
    _transform.scale(1. / _s, 1. / _s, 1.);
    _transform.translate(_Nx / 2., _Ny / 2., 0);

    // determine the solid angle subtended by a pixel in the center of the viewport
    // and use this as an approximate, representative solid angle for all pixels
    double alpha = 2. * atan(0.5 * _s / _Fe);
    double omega = alpha * alpha;

    // configure flux recorder
    instrumentFluxRecorder()->includeSurfaceBrightnessForLocal(_Nx, _Ny, omega, _s, _s, 0., 0., "length");
}

////////////////////////////////////////////////////////////////////

bool PerspectiveInstrument::hasSameSightLine(const Instrument* other) const
{
    auto perspective = dynamic_cast<const PerspectiveInstrument*>(other);
    return perspective && width() == perspective->width() && viewX() == perspective->viewX()
           && viewY() == perspective->viewY() && viewZ() == perspective->viewZ() && crossX() == perspective->crossX()
           && crossY() == perspective->crossY() && crossZ() == perspective->crossZ() && upX() == perspective->upX()
           && upY() == perspective->upY() && upZ() == perspective->upZ() && focal() == perspective->focal();
}

////////////////////////////////////////////////////////////////////

Direction PerspectiveInstrument::bfkobs(Position bfr) const
{
    // distance from launch to eye
    double Px, Py, Pz;
    bfr.cartesian(Px, Py, Pz);
    double D = norm(_Ex - Px, _Ey - Py, _Ez - Pz);

    // if the distance is very small, return something silly - the photon packet is behind the viewport anyway
    if (D < _s / 10.) return Direction();

    // otherwise return a unit vector in the direction from launch to eye
    return Direction((_Ex - Px) / D, (_Ey - Py) / D, (_Ez - Pz) / D, false);
}

////////////////////////////////////////////////////////////////////

Direction PerspectiveInstrument::bfky(Position /*bfr*/) const
{
    return _bfky;
}

////////////////////////////////////////////////////////////////////

Instrument::Detection PerspectiveInstrument::locate(Position bfr) const
{
    // get the position
    double x, y, z;
    bfr.cartesian(x, y, z);

    // transform from world to pixel coordinates
    double xp, yp, zp, wp;
    _transform.transform(x, y, z, 1., xp, yp, zp, wp);
    int i = static_cast<int>(xp / wp);
    int j = static_cast<int>(yp / wp);

    // ignore photon packets arriving outside the viewport, originating from behind the viewport,
    // or originating from very close to the viewport
    Detection detection;
    if (i >= 0 && i < _Nx && j >= 0 && j < _Ny && zp > _s / 10.)
    {
        // detect the photon packet in the appropriate pixel of the data cube and at the appropriate distance
        detection.pixel = i + _Nx * j;
        detection.distance = zp;
    }
    return detection;
}

////////////////////////////////////////////////////////////////////
