/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PeelOffInstrument.hpp"
#include "FluxRecorder.hpp"

////////////////////////////////////////////////////////////////////

void PeelOffInstrument::setupSelfBefore()
{
    Instrument::setupSelfBefore();

    instrumentFluxRecorder()->setRecordPolarization(_recordPolarization);
}

////////////////////////////////////////////////////////////////////

void PeelOffInstrument::detect(const PhotonPacket* pp, Detection detection, double extinction)
{
    instrumentFluxRecorder()->detect(pp, detection.pixel, detection.distance, extinction);
}

////////////////////////////////////////////////////////////////////
