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

bool GradientCriterion::usesNoise() const
{
    return significance() > 0.;
}

////////////////////////////////////////////////////////////////////

double GradientCriterion::excess(int m, const CellWindow& window, const CellWindow* noiseWindow,
                                 const vector<int>& neighbors) const
{
    // a cell below the floor is never flagged (a NaN value is not flagged either), nor is a cell with a value that
    // is not positive when comparing logarithms
    double value = window.mean(m);
    if (!(value >= minValue())) return 0.;
    if (logarithmic() && !(value > 0.)) return 0.;
    double field = logarithmic() ? std::log(value) : value;

    // the variance of the window-averaged field value of a cell on the scale of the comparison, i.e. the mean of the
    // squared noise of the iterations in the window divided by their number, and relative for logarithms
    bool test = noiseWindow && significance() > 0.;
    int numIterations = noiseWindow ? noiseWindow->numIterations() : 1;
    auto variance = [this, noiseWindow, numIterations](int n, double v) {
        double var = noiseWindow->mean(n) / numIterations;
        return logarithmic() ? var / (v * v) : var;
    };
    double varianceCell = test ? variance(m, value) : 0.;

    // the measure is the largest absolute difference with a face neighbor, reduced by the significance factor times
    // its standard error if requested, ignoring differences that are NaN, and neighbors with a value that is not
    // positive when comparing logarithms
    double measure = 0.;
    for (int n : neighbors)
    {
        double other = window.mean(n);
        if (logarithmic() && !(other > 0.)) continue;
        double difference = abs(field - (logarithmic() ? std::log(other) : other));
        if (test) difference -= significance() * sqrt(varianceCell + variance(n, other));
        if (difference > measure) measure = difference;
    }
    return measure > maxChange() ? measure / maxChange() : 0.;
}

////////////////////////////////////////////////////////////////////
