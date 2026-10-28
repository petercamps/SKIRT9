/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef REFINEMENTCRITERION_HPP
#define REFINEMENTCRITERION_HPP

#include "IterationHistory.hpp"
#include "RefinementCondition.hpp"

////////////////////////////////////////////////////////////////////

/** RefinementCriterion is the abstract base class for the criteria that decide which spatial
    cells are subdivided by dynamic grid refinement (see the DynamicRefinementOptions class).
    Several criteria can be configured at once, for example one per ion; a cell is subdivided if
    any of them asks for it.

    A criterion consists of three parts:
    - a \em test on a per-cell quantity, called the field of the criterion, which flags the cells
      that are too coarse; the test is implemented by the abstract subclasses GradientCriterion
      (the field differs strongly between the cell and a neighboring cell) and ValueCriterion (the
      field lies in a given range), and the field by their concrete subclasses;
    - optional \em conditions, which restrict the criterion to the cells whose state lies in given
      ranges (see the RefinementCondition class);
    - \em limits, which tell how far the criterion may subdivide a cell.

    A cell asks for subdivision, and is then called a candidate, if the criterion flags it while
    all conditions of the criterion hold and the cell has not reached a limit of the criterion. The
    test and the conditions are evaluated on the values of the field and of the condition
    quantities averaged over the iterations of the averaging window of the refinement.

    <b>Limits</b>

    The level limit is set by the \em maxExtraLevels property: the criterion does not flag a cell
    that has gained this number of levels or more beyond its level in the initial grid, whichever
    criteria caused the earlier subdivisions. For a binary tree, which splits a cell along one axis
    per level, three levels correspond to one octtree level.

    The column density limit is set by the \em minColumnDensity property: the criterion does not
    flag a cell with a column density below this value. The column density of a cell is its number
    density multiplied by its size, which is the largest extent of the cell along the coordinate
    axes. The number density is the one stored in the medium state for the medium component that
    provides the field of the criterion, or the sum over the relevant medium components, as
    described for each criterion. Because child cells inherit their parent's number density, each
    octtree subdivision halves the column density of a cell, so that a column density limit
    produces cells with a size inversely proportional to the local density.

    For both properties, a value of zero means that there is no limit. Another criterion with wider
    limits can still flag a cell that has reached the limits of this criterion.

    <b>Priority</b>

    When the cap on the number of cells does not leave room to subdivide all candidates, the
    candidates are ranked by priority. For this criterion, the priority of a candidate is the
    product of the excess, which is the ratio of the measure to the threshold for a gradient
    criterion and one for a value criterion, and the ratio of the cell's column density to the
    column density limit, or one if there is no such limit.

    <b>Responsibilities</b>

    A criterion tells whether its field changes in a given iteration loop, provides the current
    value of its field and the number density for its column density limit in a given cell, and,
    given the field averaged over the window, tells whether it flags a given cell and with which
    excess. It keeps no per-cell data of its own: the medium system accumulates its field values in
    the iteration history, and performs the refinement step. */
class RefinementCriterion : public SimulationItem
{
    ITEM_ABSTRACT(RefinementCriterion, SimulationItem, "a dynamic grid refinement criterion")

        PROPERTY_ITEM_LIST(conditions, RefinementCondition, "the conditions that restrict the criterion")
        ATTRIBUTE_REQUIRED_IF(conditions, "false")

        PROPERTY_INT(maxExtraLevels, "the maximum number of levels a cell may gain beyond its initial level "
                                     "through this criterion (zero means no limit)")
        ATTRIBUTE_MIN_VALUE(maxExtraLevels, "0")
        ATTRIBUTE_MAX_VALUE(maxExtraLevels, "99")
        ATTRIBUTE_DEFAULT_VALUE(maxExtraLevels, "0")

        PROPERTY_DOUBLE(minColumnDensity, "the column density below which this criterion does not subdivide a cell "
                                          "(zero means no limit)")
        ATTRIBUTE_QUANTITY(minColumnDensity, "numbersurfacedensity")
        ATTRIBUTE_MIN_VALUE(minColumnDensity, "[0")
        ATTRIBUTE_DEFAULT_VALUE(minColumnDensity, "0")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns true if the field of the criterion changes in the specified iteration
        loop, so that the criterion takes part in the refinement during that loop, and false
        otherwise. */
    virtual bool changesIn(IterationHistory::Loop loop) const = 0;

    /** This function returns the current value of the field of the criterion in the spatial cell
        with index \f$m\f$. The function may be called from multiple threads at the same time. */
    virtual double value(int m) const = 0;

    /** This function returns true if calculating the field value of a cell is expensive, so that
        the medium system calculates the values for all cells in parallel before accumulating
        them, and false otherwise. The default implementation returns false. */
    virtual bool hasCostlyValue() const;

    /** This function returns the number density used for the column density limit of the
        criterion in the spatial cell with index \f$m\f$. */
    virtual double numberDensity(int m) const = 0;

    /** This function returns a description of the field of the criterion, for use in log messages
        and in the description of the corresponding cell window in the iteration history. */
    virtual string fieldDescription() const = 0;

    /** This function returns true if the test of the criterion needs the face neighbors of a cell
        (see the excess() function), and false otherwise. */
    virtual bool needsNeighbors() const = 0;

    /** This function returns the excess with which the test of the criterion flags the spatial cell
        with index \f$m\f$, or zero if it does not flag the cell. The test is applied to the field
        values averaged over the window, which the function obtains from the specified cell window
        for any cell index. If needsNeighbors() returns true, the vector \em neighbors holds the
        indices of the face neighbors of the cell; otherwise it is empty. The function does not
        consider the conditions or the limits of the criterion. It may be called from multiple
        threads at the same time. */
    virtual double excess(int m, const CellWindow& window, const vector<int>& neighbors) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
