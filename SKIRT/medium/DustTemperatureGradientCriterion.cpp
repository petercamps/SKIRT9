/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustTemperatureGradientCriterion.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"

////////////////////////////////////////////////////////////////////

void DustTemperatureGradientCriterion::setupSelfAfter()
{
    GradientCriterion::setupSelfAfter();

    // the medium system is not set up, because the criterion is part of it; the lists of medium indices are available
    _ms = find<MediumSystem>(false);
    if (_ms->dustMediumIndices().empty())
        throw FATALERROR("A dust temperature refinement criterion requires a dust medium component");
    if (!find<Configuration>()->hasRadiationField())
        throw FATALERROR("A dust temperature refinement criterion requires the radiation field to be stored");
}

////////////////////////////////////////////////////////////////////

bool DustTemperatureGradientCriterion::changesIn(IterationHistory::Loop loop) const
{
    return loop != IterationHistory::Loop::None;
}

////////////////////////////////////////////////////////////////////

double DustTemperatureGradientCriterion::value(int m) const
{
    return _ms->indicativeDustTemperature(m);
}

////////////////////////////////////////////////////////////////////

bool DustTemperatureGradientCriterion::hasCostlyValue() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

double DustTemperatureGradientCriterion::noise(int m) const
{
    // the change of the temperature when the radiation field changes by its relative noise; a cell without any
    // contributions to the radiation field has infinite noise, and a noise that is not tracked is zero
    double sigma = _ms->radiationFieldNoise(m);
    if (!(sigma > 0.) || std::isinf(sigma)) return sigma;
    return abs(_ms->indicativeDustTemperature(m, 1. + sigma) - _ms->indicativeDustTemperature(m));
}

////////////////////////////////////////////////////////////////////

double DustTemperatureGradientCriterion::numberDensity(int m) const
{
    double n = 0.;
    for (int h : _ms->dustMediumIndices()) n += _ms->numberDensity(m, h);
    return n;
}

////////////////////////////////////////////////////////////////////

string DustTemperatureGradientCriterion::fieldDescription() const
{
    return "indicative dust temperature";
}

////////////////////////////////////////////////////////////////////
