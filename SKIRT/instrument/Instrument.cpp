/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Instrument.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "FluxRecorder.hpp"

////////////////////////////////////////////////////////////////////

void Instrument::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // select "local" or default wavelength grid
    auto config = find<Configuration>();
    _instrumentWavelengthGrid = config->wavelengthGrid(wavelengthGrid());

    // discover details about the simulation
    bool hasMedium = config->hasMedium();
    bool hasMediumEmission = config->hasSecondaryEmission() || config->scatteringEmulatesSecondaryEmission();

    // partially configure the flux recorder
    _recorder = new FluxRecorder(this);
    _recorder->setSimulationInfo(instrumentName(), hasMedium, hasMediumEmission);
    _recorder->setWavelengthGrid(instrumentWavelengthGrid());
    _recorder->setUserFlags(_recordComponents, _numScatteringLevels, _recordStatistics);
}

////////////////////////////////////////////////////////////////////

void Instrument::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // finalize configuration of the flux recorder
    _recorder->finalizeConfiguration();
}

////////////////////////////////////////////////////////////////////

Instrument::~Instrument()
{
    delete _recorder;
}

////////////////////////////////////////////////////////////////////

void Instrument::configureRecorderDistance(double distance)
{
    if (distance > 0.)
    {
        _recorder->setRestFrameDistance(distance);
    }
    else
    {
        auto config = find<Configuration>();
        if (config->redshift() > 0.)
        {
            _recorder->setObserverFrameRedshift(config->redshift(), config->angularDiameterDistance(),
                                                config->luminosityDistance());
        }
        else
        {
            throw FATALERROR("Instrument distance and model redshift are both zero");
        }
    }
}

////////////////////////////////////////////////////////////////////

std::string Instrument::itemName() const
{
    return instrumentName();
}

////////////////////////////////////////////////////////////////////

void Instrument::flush()
{
    _recorder->flush();
}

////////////////////////////////////////////////////////////////////

void Instrument::write()
{
    _recorder->calibrateAndWrite();
}

////////////////////////////////////////////////////////////////////

bool Instrument::recordsWavelength(double lambda) const
{
    return _recorder->recordsWavelength(lambda);
}

////////////////////////////////////////////////////////////////////
