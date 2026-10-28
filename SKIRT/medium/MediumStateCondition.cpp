/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumStateCondition.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void MediumStateCondition::setupSelfAfter()
{
    RefinementCondition::setupSelfAfter();

    if (maxValue() < minValue()) throw FATALERROR("The maximum value of a condition cannot be below its minimum value");
    _selection.select(this, variable(), false);
}

////////////////////////////////////////////////////////////////////

double MediumStateCondition::value(int m) const
{
    return _selection.value(m);
}

////////////////////////////////////////////////////////////////////

bool MediumStateCondition::holds(double value) const
{
    return value >= minValue() && value <= maxValue();
}

////////////////////////////////////////////////////////////////////

string MediumStateCondition::description() const
{
    return "condition variable " + _selection.description();
}

////////////////////////////////////////////////////////////////////
