/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HISTORYSCOPE_HPP
#define HISTORYSCOPE_HPP

#include "IterationHistory.hpp"

////////////////////////////////////////////////////////////////////

/** A HistoryScope gives a client access to the series kept in the IterationHistory on behalf of a
    given simulation item, without the client knowing that item. The medium system passes a scope
    bound to a medium component to the convergence function of the material mix of that component,
    so that the series of the mix are keyed on the medium component rather than on the mix. As a
    result, the mix need not know its component, and its keys remain stable when the component uses
    a material mix family. A scope is a lightweight handle that can be copied freely. */
class HistoryScope
{
public:
    /** This constructor creates a scope for the specified history and item. */
    HistoryScope(IterationHistory* history, const SimulationItem* item) : _history{history}, _item{item} {}

    /** This function returns the scalar series with the key (item, id), as described for
        IterationHistory::scalarSeries(). */
    ScalarSeries& scalarSeries(int id, int depth, IterationHistory::Lifetime lifetime, string description,
                               string quantity = string()) const
    {
        return _history->scalarSeries(_item, id, depth, lifetime, description, quantity);
    }

    /** This function returns the existing scalar series with the key (item, id), as described for
        IterationHistory::series(). A material mix uses this function to retrieve an aggregate
        series declared on its behalf by the medium system. */
    const ScalarSeries& series(int id) const { return _history->series(_item, id); }

    /** This function returns the current iteration loop, as described for
        IterationHistory::loop(). */
    IterationHistory::Loop loop() const { return _history->loop(); }

    /** This function returns the one-based index of the current iteration, as described for
        IterationHistory::iteration(). */
    int iteration() const { return _history->iteration(); }

private:
    IterationHistory* _history{nullptr};
    const SimulationItem* _item{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
