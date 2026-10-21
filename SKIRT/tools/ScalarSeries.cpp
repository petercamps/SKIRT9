/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ScalarSeries.hpp"

////////////////////////////////////////////////////////////////////

ScalarSeries::ScalarSeries(SeriesKey key, int depth, string description, string quantity)
    : _key{key}, _description{description}, _quantity{quantity}, _values(max(depth, 1))
{}

////////////////////////////////////////////////////////////////////

void ScalarSeries::set(double value)
{
    _values[0] = value;
}

////////////////////////////////////////////////////////////////////

bool ScalarSeries::has(int lag) const
{
    return lag >= 0 && lag < depth() && _values[lag].has_value();
}

////////////////////////////////////////////////////////////////////

double ScalarSeries::value(int lag) const
{
    return has(lag) ? *_values[lag] : std::numeric_limits<double>::quiet_NaN();
}

////////////////////////////////////////////////////////////////////

double ScalarSeries::relativeChange(int lag) const
{
    if (!has(0) || !has(lag)) return std::numeric_limits<double>::quiet_NaN();
    double current = value(0);
    double previous = value(lag);
    if (current == previous) return 0.;
    if (current == 0.) return std::numeric_limits<double>::infinity();
    return abs(current - previous) / abs(current);
}

////////////////////////////////////////////////////////////////////

bool ScalarSeries::isStable(int count, double tolerance) const
{
    if (count < 1 || count > depth()) return false;
    for (int lag = 0; lag != count; ++lag)
        if (!has(lag)) return false;
    for (int lag = 0; lag + 1 < count; ++lag)
        if (abs(value(lag) - value(lag + 1)) > tolerance) return false;
    return true;
}

////////////////////////////////////////////////////////////////////

void ScalarSeries::grow(int newDepth)
{
    if (newDepth > depth()) _values.resize(newDepth);
}

////////////////////////////////////////////////////////////////////

void ScalarSeries::setAggregationRule(const AggregationRule& rule)
{
    _hasRule = true;
    _rule = rule;
}

////////////////////////////////////////////////////////////////////

void ScalarSeries::advance()
{
    for (int lag = depth() - 1; lag > 0; --lag) _values[lag] = _values[lag - 1];
    _values[0].reset();
}

////////////////////////////////////////////////////////////////////

void ScalarSeries::clear()
{
    for (auto& value : _values) value.reset();
}

////////////////////////////////////////////////////////////////////
