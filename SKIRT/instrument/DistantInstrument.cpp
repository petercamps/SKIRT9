/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DistantInstrument.hpp"
#include "FluxRecorder.hpp"

////////////////////////////////////////////////////////////////////

void DistantInstrument::setupSelfBefore()
{
    PeelOffInstrument::setupSelfBefore();

    // pass the observer angles to the flux recorder
    instrumentFluxRecorder()->setObserverAngles(_inclination, _azimuth, _roll);

    // configure the flux recorder with the appropriate frame and distances
    configureRecorderDistance(distance());

    // calculate sine and cosine for our angles
    double costheta = cos(_inclination);
    double sintheta = sin(_inclination);
    double cosphi = cos(_azimuth);
    double sinphi = sin(_azimuth);
    double cosomega = cos(_roll);
    double sinomega = sin(_roll);

    // calculate relevant directions
    _bfkobs = Direction(_inclination, _azimuth);
    _bfky = Direction(-cosphi * costheta * cosomega - sinphi * sinomega,
                      -sinphi * costheta * cosomega + cosphi * sinomega, +sintheta * cosomega, false);
}

////////////////////////////////////////////////////////////////////

bool DistantInstrument::hasSameSightLine(const PeelOffInstrument* other) const
{
    // the roll angle determines the orientation of the reference frame for the polarization state; it
    // does not affect the intensity and thus matters only if any of the instruments records polarization
    auto distant = dynamic_cast<const DistantInstrument*>(other);
    return distant && inclination() == distant->inclination() && azimuth() == distant->azimuth()
           && (roll() == distant->roll() || (!recordPolarization() && !distant->recordPolarization()));
}

////////////////////////////////////////////////////////////////////

Direction DistantInstrument::bfkobs(Position /*bfr*/) const
{
    return _bfkobs;
}

////////////////////////////////////////////////////////////////////

Direction DistantInstrument::bfky(Position /*bfr*/) const
{
    return _bfky;
}

////////////////////////////////////////////////////////////////////
