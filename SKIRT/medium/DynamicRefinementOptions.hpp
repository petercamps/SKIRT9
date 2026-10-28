/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DYNAMICREFINEMENTOPTIONS_HPP
#define DYNAMICREFINEMENTOPTIONS_HPP

#include "RefinementCriterion.hpp"

////////////////////////////////////////////////////////////////////

/** The DynamicRefinementOptions class holds the settings for dynamic grid refinement, i.e. the
    subdivision of spatial cells between the iterations of a simulation that iterates over primary
    or secondary emission. It is configured through the \em dynamicRefinementOptions property of a
    tree spatial grid (see the TreeSpatialGrid class), which is the only type of grid that supports
    dynamic refinement.

    <b>Overview</b>

    Some simulations depend on spatial structure that is not known when the spatial grid is built.
    For example, an ionization front, or the thin layer in which a given ion is concentrated, only
    emerges once the radiation field has been calculated. Dynamic refinement subdivides cells
    between iterations, wherever the medium state or the radiation field calculated so far shows
    that a cell is too coarse: because a quantity varies too steeply across the cell, or because
    the cell holds material in a state that must be resolved down to a given scale. The simulation
    continues on the finer grid, and the iterations proceed until both the medium state and the
    grid have stopped changing. Cells are only ever subdivided, never merged.

    A simulation with dynamic refinement proceeds as follows:
    - The tree is constructed as usual by the configured tree policies. This initial grid should
      already be adequate for the density distribution and, where applicable, resolve the expected
      areas of special interest.
    - The iterations start. After the radiation field and the medium state have been updated at the
      end of each iteration, the refinement step records, for each configured refinement criterion,
      the value of its field in each cell.
    - At the end of each averaging window of a few iterations, each refinement criterion evaluates
      the window-averaged fields. The cells that a criterion flags as too coarse, and that have not
      yet reached the limits of that criterion, are subdivided (see the RefinementCriterion class).
    - The iterations continue on the refined grid. The loop terminates once the medium state has
      converged and the refinement has settled, meaning that an averaging window has passed without
      leading to further subdivision.

    <b>Schedule</b>

    Dynamic refinement can take place in each of the iteration loops: primary emission iterations,
    secondary emission iterations, and merged primary and secondary emission iterations. In a given
    loop, only the criteria whose field changes in that loop take part. If no criterion takes part,
    the refinement step does nothing in that loop. In each loop, the refinement starts afresh and
    goes through the following phases:
    - During the first \em numIterationsBeforeRefinement iterations, refinement is not considered
      at all, so that the first decision is based on a settled radiation field.
    - The refinement fields are then accumulated for \em numAveragedIterations iterations, after
      which the criteria evaluate the window-averaged values, and the refinement step decides which
      cells to subdivide. Averaging over several iterations reduces the Monte Carlo noise that
      reaches the decision.
    - If cells are subdivided, which is called a refinement round, the next \em
      numIterationsAfterRefinement iterations are left out of the averaging, because the field is
      still adjusting to the changed grid. A new window starts after that. If no cells are
      subdivided, a new window starts right away.

    A decision does not lead to a refinement round, even though cells ask for subdivision, if the
    loop has already performed \em maxRounds refinement rounds, if the \em maxCells cap leaves no
    room for subdividing even a single cell, or if fewer than \em numIterationsAfterRefinement +
    \em numAveragedIterations iterations remain before the loop reaches its maximum number of
    iterations. In each of these cases, a warning reports the number of candidates that were left
    unsplit.

    The iteration loop is considered converged only if the medium state has converged and the
    refinement has settled. The refinement has settled in an iteration if a decision was taken in
    that iteration and did not lead to a refinement round. As a result, each refinement round costs
    at least \em numIterationsAfterRefinement + \em numAveragedIterations extra iterations. If the
    loop ends because it reaches the maximum number of iterations while the refinement has not
    settled, a warning is issued.

    <b>Cell cap</b>

    The \em maxCells property sets a hard cap on the total number of cells, or zero for no cap.
    When the cap prevents subdividing all candidates, the candidates with the highest priority are
    subdivided first (see the RefinementCriterion class), and a warning reports how many cells were
    left unsplit. Among candidates with equal priority, the larger cells go first, and among those
    the cells with the lower cell index. If a cap is configured, the memory for that number of
    cells is reserved at the end of setup, so that the memory need of the run is fixed from the
    start. Setup logs the amount of memory reserved for each process, and issues a warning if the
    memory of the process with all of that capacity in use would exceed the physical memory of the
    host. It also warns if the number of photon packets per cell launched in an iteration of a loop
    in which refinement takes place would be below 100 for a grid with the maximum number of cells.

    <b>Initial state of child cells</b>

    A child cell inherits its parent's medium state and material mix, so that the total mass of
    each medium component is conserved, and the next iteration starts from the previous solution.
    Intensive and categorical state variables are copied to each child, and extensive state
    variables are divided among the children in proportion to their volume. The radiation field
    and the luminosity absorbed by the dust in the parent are divided among the children in
    proportion to their volume as well, so that each child starts with its parent's mean intensity.
    Inheriting the state means that the refinement does not resolve any density structure within
    the parent cell. */
class DynamicRefinementOptions : public SimulationItem
{
    ITEM_CONCRETE(DynamicRefinementOptions, SimulationItem, "a set of options for dynamic grid refinement")

        PROPERTY_INT(maxCells, "the maximum number of cells in the grid (zero means no maximum)")
        ATTRIBUTE_MIN_VALUE(maxCells, "0")
        ATTRIBUTE_MAX_VALUE(maxCells, "2000000000")
        ATTRIBUTE_DEFAULT_VALUE(maxCells, "0")

        PROPERTY_INT(numIterationsBeforeRefinement,
                     "the number of initial iterations of each loop during which refinement is not considered")
        ATTRIBUTE_MIN_VALUE(numIterationsBeforeRefinement, "0")
        ATTRIBUTE_MAX_VALUE(numIterationsBeforeRefinement, "1000")
        ATTRIBUTE_DEFAULT_VALUE(numIterationsBeforeRefinement, "5")

        PROPERTY_INT(numAveragedIterations,
                     "the number of iterations over which the refinement fields are averaged before each decision")
        ATTRIBUTE_MIN_VALUE(numAveragedIterations, "1")
        ATTRIBUTE_MAX_VALUE(numAveragedIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(numAveragedIterations, "3")

        PROPERTY_INT(numIterationsAfterRefinement,
                     "the number of iterations after each refinement round left out of the averaging")
        ATTRIBUTE_MIN_VALUE(numIterationsAfterRefinement, "0")
        ATTRIBUTE_MAX_VALUE(numIterationsAfterRefinement, "1000")
        ATTRIBUTE_DEFAULT_VALUE(numIterationsAfterRefinement, "1")
        ATTRIBUTE_DISPLAYED_IF(numIterationsAfterRefinement, "Level3")

        PROPERTY_INT(maxRounds, "the maximum number of refinement rounds in each loop (zero means no maximum)")
        ATTRIBUTE_MIN_VALUE(maxRounds, "0")
        ATTRIBUTE_MAX_VALUE(maxRounds, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxRounds, "0")

        PROPERTY_ITEM_LIST(criteria, RefinementCriterion, "the refinement criteria")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
