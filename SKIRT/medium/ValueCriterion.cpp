/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ValueCriterion.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void ValueCriterion::setupSelfBefore()
{
    RefinementCriterion::setupSelfBefore();

    if (maxValue() < minValue())
        throw FATALERROR("The maximum value of a refinement criterion cannot be below its minimum value");
    if (!maxExtraLevels() && !minColumnDensity())
        throw FATALERROR("A value refinement criterion requires a level limit or a column density limit");
}

////////////////////////////////////////////////////////////////////

bool ValueCriterion::needsNeighbors() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

double ValueCriterion::excess(int m, const CellWindow& window, const vector<int>& /*neighbors*/) const
{
    double value = window.mean(m);
    return value >= minValue() && value <= maxValue() ? 1. : 0.;
}

////////////////////////////////////////////////////////////////////
