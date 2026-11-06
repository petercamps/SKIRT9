/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumStateGradientCriterion.hpp"
#include "DynamicStateOptions.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"

////////////////////////////////////////////////////////////////////

void MediumStateGradientCriterion::setupSelfAfter()
{
    GradientCriterion::setupSelfAfter();
    _selection.select(this, variable(), true);

    // the noise test needs a recipe that can recalculate the variable for a scaled radiation field
    if (significance() > 0.)
    {
        // the medium system is not set up, because the criterion is part of it; its media have been set up, however
        _ms = find<MediumSystem>(false);
        if (_ms->dynamicStateOptions())
        {
            for (auto recipe : _ms->dynamicStateOptions()->recipes())
            {
                if (recipe->offersScaledFieldStateVariable(_selection.mediumIndex(), _selection.customIndex()))
                {
                    _recipe = recipe;
                    break;
                }
            }
        }
        if (!_recipe)
            find<Log>()->warning("The refinement criterion on medium state variable " + variable()
                                 + " ignores its significance property, because no dynamic medium state recipe "
                                   "provides the noise of the variable");
    }
}

////////////////////////////////////////////////////////////////////

bool MediumStateGradientCriterion::changesIn(IterationHistory::Loop loop) const
{
    return _selection.changesIn(loop);
}

////////////////////////////////////////////////////////////////////

double MediumStateGradientCriterion::value(int m) const
{
    return _selection.value(m);
}

////////////////////////////////////////////////////////////////////

double MediumStateGradientCriterion::numberDensity(int m) const
{
    return _selection.numberDensity(m);
}

////////////////////////////////////////////////////////////////////

bool MediumStateGradientCriterion::usesNoise() const
{
    return _recipe != nullptr;
}

////////////////////////////////////////////////////////////////////

bool MediumStateGradientCriterion::hasCostlyValue() const
{
    return usesNoise();
}

////////////////////////////////////////////////////////////////////

double MediumStateGradientCriterion::noise(int m) const
{
    // a cell without any contributions to the radiation field has infinite noise, and a noise that is not tracked
    // is zero
    double sigma = _ms->radiationFieldNoise(m);
    if (!(sigma > 0.) || std::isinf(sigma)) return sigma;

    // the difference between the values of the variable for the radiation field scaled up and down by k standard
    // errors, divided by 2k
    double k = significance();
    int h = _selection.mediumIndex();
    int i = _selection.customIndex();
    double lo = _recipe->scaledFieldStateVariable(m, h, i, max(0., 1. - k * sigma));
    double hi = _recipe->scaledFieldStateVariable(m, h, i, 1. + k * sigma);
    return abs(hi - lo) / (2. * k);
}

////////////////////////////////////////////////////////////////////

string MediumStateGradientCriterion::fieldDescription() const
{
    return _selection.description();
}

////////////////////////////////////////////////////////////////////
