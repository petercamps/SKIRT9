/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ApertureInstrument.hpp"

////////////////////////////////////////////////////////////////////

void ApertureInstrument::setupSelfBefore()
{
    DistantInstrument::setupSelfBefore();

    // precalculate information needed by locate() function
    _radius2 = radius() * radius();
    _costheta = cos(inclination());
    _sintheta = sin(inclination());
    _cosphi = cos(azimuth());
    _sinphi = sin(azimuth());
}

////////////////////////////////////////////////////////////////////

Instrument::Detection ApertureInstrument::locate(Position bfr) const
{
    Detection detection;

    // if the instrument has an aperture
    if (_radius2)
    {
        // get the position
        double x, y, z;
        bfr.cartesian(x, y, z);

        // transform to detector coordinates using inclination and azimuth
        // but without performing the roll, which would not alter the radius
        double xpp = -_sinphi * x + _cosphi * y;
        double ypp = -_cosphi * _costheta * x - _sinphi * _costheta * y + _sintheta * z;

        // if the position projects outside of the aperture, the photon packet is not detected
        double radius2 = xpp * xpp + ypp * ypp;
        if (radius2 > _radius2) return detection;
    }

    // otherwise, the photon packet is detected in the single "pixel"
    detection.pixel = 0;
    return detection;
}

////////////////////////////////////////////////////////////////////
