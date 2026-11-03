/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RadiationFieldGradientCriterion.hpp"
#include "Configuration.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"

////////////////////////////////////////////////////////////////////

void RadiationFieldGradientCriterion::setupSelfAfter()
{
    GradientCriterion::setupSelfAfter();

    // the medium system is not set up, because the criterion is part of it
    _ms = find<MediumSystem>(false);
    auto config = find<Configuration>();
    if (!config->hasRadiationField())
        throw FATALERROR("A radiation field refinement criterion requires the radiation field to be stored");

    // determine the bins of the radiation field wavelength grid in the integration range
    _wavelengthGrid = config->radiationFieldWLG();
    int numBins = _wavelengthGrid->numBins();
    for (int ell = 0; ell != numBins; ++ell)
    {
        double lambda = _wavelengthGrid->wavelength(ell);
        if (lambda >= minWavelength() && lambda <= maxWavelength()) _bins.push_back(ell);
    }
    if (_bins.empty())
        throw FATALERROR("The integration range of a radiation field refinement criterion contains no bin of the "
                         "radiation field wavelength grid");
}

////////////////////////////////////////////////////////////////////

bool RadiationFieldGradientCriterion::changesIn(IterationHistory::Loop loop) const
{
    return loop != IterationHistory::Loop::None;
}

////////////////////////////////////////////////////////////////////

double RadiationFieldGradientCriterion::value(int m) const
{
    Array Jv = _ms->meanIntensity(m);
    double J = 0.;
    for (int ell : _bins) J += Jv[ell] * _wavelengthGrid->effectiveWidth(ell);
    return J;
}

////////////////////////////////////////////////////////////////////

double RadiationFieldGradientCriterion::noise(int m) const
{
    // a cell without any contributions to the radiation field has infinite noise, and a noise that is not tracked
    // is zero
    double sigma = _ms->radiationFieldNoise(m);
    if (!(sigma > 0.) || std::isinf(sigma)) return sigma;
    return value(m) * sigma;
}

////////////////////////////////////////////////////////////////////

double RadiationFieldGradientCriterion::numberDensity(int m) const
{
    double n = 0.;
    for (int h = 0; h != _ms->numMedia(); ++h) n += _ms->numberDensity(m, h);
    return n;
}

////////////////////////////////////////////////////////////////////

string RadiationFieldGradientCriterion::fieldDescription() const
{
    return "integrated mean intensity";
}

////////////////////////////////////////////////////////////////////
