/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GradientCriterion.hpp"

////////////////////////////////////////////////////////////////////

bool GradientCriterion::needsNeighbors() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

double GradientCriterion::excess(int m, const CellWindow& window, const vector<int>& neighbors) const
{
    // a cell below the floor is never flagged (a NaN value is not flagged either)
    double value = window.mean(m);
    if (!(value >= minValue())) return 0.;

    // the measure is the largest absolute difference with a face neighbor (NaN differences are ignored)
    double measure = 0.;
    for (int n : neighbors)
    {
        double difference = abs(value - window.mean(n));
        if (difference > measure) measure = difference;
    }
    return measure > maxChange() ? measure / maxChange() : 0.;
}

////////////////////////////////////////////////////////////////////
