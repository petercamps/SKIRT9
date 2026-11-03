/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumStateGradientCriterion.hpp"
#include "Log.hpp"

////////////////////////////////////////////////////////////////////

void MediumStateGradientCriterion::setupSelfAfter()
{
    GradientCriterion::setupSelfAfter();
    _selection.select(this, variable(), true);

    // the noise of a medium state variable is not known, so the noise test cannot be performed
    if (significance() > 0.)
        find<Log>()->warning("The refinement criterion on medium state variable " + variable()
                             + " ignores its significance property, because the noise of the variable is not known");
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
    return false;
}

////////////////////////////////////////////////////////////////////

string MediumStateGradientCriterion::fieldDescription() const
{
    return _selection.description();
}

////////////////////////////////////////////////////////////////////
