/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SCALARSERIES_HPP
#define SCALARSERIES_HPP

#include "AggregationRule.hpp"
#include "SeriesKey.hpp"
#include <optional>

////////////////////////////////////////////////////////////////////

/** A ScalarSeries holds the values of a single global quantity over the most recent iterations,
    up to its depth, on behalf of a client of the IterationHistory. All scalar series in the
    history advance together, one slot per iteration, so that "one iteration back" means the same
    thing for every series. During an iteration, the client sets the value for the current slot; a
    slot that is not set in a given iteration remains empty, and the functions that compare values
    report a missing value rather than using a stale one.

    A series carries a human-readable description and, optionally, a physical quantity, as a
    quantity name known to the Units class, so that its values can be logged and probed. An
    aggregate series also carries an aggregation rule, according to which the medium system sets
    its value. Instances are created only by the IterationHistory class. */
class ScalarSeries
{
    //============= Construction =============

private:
    // only the IterationHistory class creates and manages series
    friend class IterationHistory;

    /** This constructor creates an empty series with the specified key, depth (at least one),
        description, and quantity. */
    ScalarSeries(SeriesKey key, int depth, string description, string quantity);

    //============= Setting =============

public:
    /** This function sets the value for the current iteration, replacing any value set earlier in
        the same iteration. The value must be identical on all processes. */
    void set(double value);

    //============= Querying =============

public:
    /** This function returns true if a value was set \em lag iterations back, where zero indicates
        the current iteration. It returns false if no value was set in that iteration, if the
        series was cleared since then, or if \em lag is negative or not smaller than the depth of
        the series. */
    bool has(int lag = 0) const;

    /** This function returns the value set \em lag iterations back, or NaN if has(lag) returns
        false. */
    double value(int lag = 0) const;

    /** This function returns the change of the current value relative to the value \em lag
        iterations back, i.e. \f$|v_0 - v_\mathrm{lag}| / |v_0|\f$. It returns zero if both values
        are equal, infinity if \f$v_0\f$ is zero and \f$v_\mathrm{lag}\f$ is not, and NaN if either
        value is missing. */
    double relativeChange(int lag = 1) const;

    /** This function returns true if the \em count most recent values (including the current one)
        are all available and each of them differs from the next by at most \em tolerance, and
        false otherwise, including if \em count exceeds the depth of the series. */
    bool isStable(int count, double tolerance) const;

    /** This function returns the key of the series. */
    const SeriesKey& key() const { return _key; }

    /** This function returns the human-readable description of the series. */
    const string& description() const { return _description; }

    /** This function returns the depth of the series, i.e. the number of iterations for which it
        keeps a value, including the current one. */
    int depth() const { return static_cast<int>(_values.size()); }

    /** This function returns the physical quantity of the values, as a quantity name known to the
        Units class, or the empty string if the values are dimensionless. */
    const string& quantity() const { return _quantity; }

    /** This function returns a pointer to the aggregation rule if this is an aggregate series, or
        the null pointer otherwise. */
    const AggregationRule* aggregationRule() const { return _hasRule ? &_rule : nullptr; }

    //============= Managing (by the history) =============

private:
    /** This function grows the series to the specified depth, preserving its values, if the depth
        is larger than the current depth. */
    void grow(int newDepth);

    /** This function attaches the specified aggregation rule to the series. */
    void setAggregationRule(const AggregationRule& rule);

    /** This function advances the series by one slot, discarding the oldest value and leaving the
        new current slot empty. */
    void advance();

    /** This function removes all values from the series. */
    void clear();

    //============= Data members =============

private:
    SeriesKey _key;
    string _description;
    string _quantity;
    bool _hasRule{false};
    AggregationRule _rule;
    vector<std::optional<double>> _values;  // the value for each lag, empty if not set
};

////////////////////////////////////////////////////////////////////

#endif
