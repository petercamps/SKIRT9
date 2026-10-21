/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ITERATIONHISTORY_HPP
#define ITERATIONHISTORY_HPP

#include "CellWindow.hpp"
#include "ScalarSeries.hpp"
#include "SimulationItem.hpp"
#include <map>

////////////////////////////////////////////////////////////////////

/** IterationHistory holds the historical data used by the iteration loops of a simulation to
    evaluate convergence. The MonteCarloSimulation creates a single instance as a child, so that any
    simulation item can locate it with find<IterationHistory>(). Like the Configuration class, it is
    based on SimulationItem but not discoverable, because it is not intended to be configured by
    the user.

    The history holds two kinds of series, each identified by a SeriesKey: a ScalarSeries holds the
    values of a global quantity over the most recent iterations, and a CellWindow holds the running
    sum of a quantity for each spatial cell. The history knows nothing about the quantities it
    stores or about the criteria that use them. All series are created on demand: if no client
    requests a series, the history holds no data. Declaring a series is idempotent, so that a
    client can simply declare its series each time it needs it. A material mix accesses the series
    kept on behalf of its medium component through a HistoryScope.

    <b>Iteration alignment</b>

    All scalar series advance together, one slot per iteration. The iteration loops of the
    simulation call beginLoop() at the start of each loop (primary, secondary, or merged), and
    beginIteration() at the start of each iteration. During an iteration, a client sets the value
    for the current slot; the last value set counts.

    <b>Lifetime</b>

    Each series has one of three lifetimes: \em loop lifetime, cleared at the start of each
    iteration loop; \em grid lifetime, cleared at the start of each iteration loop and also at the
    start of the first iteration after the spatial grid has changed, for values that cannot be
    compared between different grids; and \em simulation lifetime, kept for the complete
    simulation. Without dynamic grid refinement, the grid never changes, and grid lifetime is
    equivalent to loop lifetime.

    <b>Aggregates</b>

    An aggregate series is a scalar series with an aggregation rule, which specifies a custom
    medium state variable of a medium component. The medium system sets the value of every
    aggregate series to the volume integral of that variable over all cells, at the end of its
    setup (recording the initial state) and after each update of the medium state. Aggregate series
    have simulation lifetime, and must be declared during setup, before the medium system records
    the initial values.

    <b>Parallelization</b>

    The history is replicated on each process. All processes declare the same series in the same
    order and record the same values, so that the histories remain identical without any
    communication. A recorded value must therefore be identical on all processes, which is the case
    for values calculated from synchronized data. The listing functions return the series in the
    order in which they were declared, which is thus the same on all processes and in all runs. */
class IterationHistory : public SimulationItem
{
    //============= Data types =============

public:
    /** This enumeration lists the iteration loops of a simulation, or \c None before the first
        loop has started. */
    enum class Loop { None, Primary, Secondary, Merged };

    /** This enumeration lists the lifetimes of a series, as described in the class header. */
    enum class Lifetime { Loop, Grid, Simulation };

    //============= Construction =============

public:
    /** This constructor creates an empty history that is hooked up as a child to the specified
        parent in the simulation hierarchy, so that it will automatically be deleted. */
    explicit IterationHistory(SimulationItem* parent);

    //======== Iteration events, called by the iteration loops ========

public:
    /** This function starts a new iteration loop: it clears all series and cell windows with loop
        or grid lifetime, and resets the iteration index to zero. */
    void beginLoop(Loop loop);

    /** This function starts a new iteration: it increments the iteration index and advances all
        scalar series by one slot, leaving the new current slot empty. If the grid has changed since
        the previous call (see appendCells()), the function first clears all series and cell
        windows with grid lifetime. */
    void beginIteration();

    /** This function returns the current iteration loop, or the last loop if it has finished, or
        \c None before the first loop has started. */
    Loop loop() const { return _loop; }

    /** This function returns the one-based index of the current iteration within the current loop,
        or zero before the first iteration of the loop. */
    int iteration() const { return _iteration; }

    //======== Declaring series ========

public:
    /** This function returns the scalar series with the key (item, id), creating it if needed. If
        the series exists with a smaller depth, it grows to the requested depth, preserving its
        values. The function throws a fatal error if the series exists with a different lifetime
        or with an aggregation rule. The optional quantity is a quantity name known to the Units
        class; it is used only for output. The returned reference remains valid for the lifetime of
        the history. */
    ScalarSeries& scalarSeries(const SimulationItem* item, int id, int depth, Lifetime lifetime, string description,
                               string quantity = string());

    /** This function returns the aggregate series with the key (item, id), creating it if needed
        with the specified aggregation rule and simulation lifetime. If the series exists with a
        smaller depth, it grows to the requested depth. The function throws a fatal error if the
        series exists without an aggregation rule or with a different one. The medium system
        records the value of the series after setup and after each medium state update, so the
        series must be declared during setup. The returned reference remains valid for the lifetime
        of the history. */
    ScalarSeries& aggregateSeries(const SimulationItem* item, int id, const AggregationRule& rule, int depth = 2,
                                  string description = string(), string quantity = string());

    /** This function returns the cell window with the key (item, id), creating it for the
        specified number of cells if needed. The function throws a fatal error if the window exists
        with a different lifetime. The returned reference remains valid for the lifetime of the
        history. */
    CellWindow& cellWindow(const SimulationItem* item, int id, int numCells, Lifetime lifetime, string description);

    //======== Retrieving series ========

public:
    /** This function returns the existing scalar series with the key (item, id), for example an
        aggregate series declared during setup. It throws a fatal error if there is no such series.
        */
    ScalarSeries& series(const SimulationItem* item, int id);

    /** This function returns all aggregate series, in the order in which they were declared, so
        that the medium system can record their values. */
    vector<ScalarSeries*> aggregates();

    //======== Dynamic grid refinement ========

public:
    /** This function grows all cell windows by one cell for each element of \em parentCells. The
        new cells receive consecutive indices following the existing cells, and start with the
        running sum of the indicated parent cell. The function also marks the grid as changed, so
        that the next call to beginIteration() clears all series and cell windows with grid
        lifetime. */
    void appendCells(const vector<int>& parentCells);

    //======== Listing ========

public:
    /** This function returns all scalar series, in the order in which they were declared, for
        logging and probing. */
    vector<const ScalarSeries*> allScalarSeries() const;

    /** This function returns all cell windows, in the order in which they were declared. */
    vector<const CellWindow*> allCellWindows() const;

    //======== Private helpers ========

private:
    /** This function clears the scalar series and cell windows with the specified lifetime. */
    void clear(Lifetime lifetime);

    //======================== Data Members ========================

private:
    // a scalar series or cell window, with its lifetime
    struct ScalarEntry
    {
        ScalarSeries series;
        Lifetime lifetime;
    };
    struct WindowEntry
    {
        CellWindow window;
        Lifetime lifetime;
    };

    Loop _loop{Loop::None};
    int _iteration{0};
    bool _gridChanged{false};  // true if cells were appended since the previous call to beginIteration()

    // the series indexed on their key (node-based containers, so that references remain valid),
    // and pointers to the same series in the order of declaration
    std::map<SeriesKey, ScalarEntry> _scalars;
    std::map<SeriesKey, WindowEntry> _windows;
    vector<ScalarEntry*> _scalarOrder;
    vector<WindowEntry*> _windowOrder;
};

////////////////////////////////////////////////////////////////////

#endif
