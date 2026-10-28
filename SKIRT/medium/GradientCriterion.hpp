/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GRADIENTCRITERION_HPP
#define GRADIENTCRITERION_HPP

#include "RefinementCriterion.hpp"

////////////////////////////////////////////////////////////////////

/** GradientCriterion is the abstract base class for dynamic grid refinement criteria that flag a
    spatial cell if its field differs strongly from the field in a neighboring cell (see the
    RefinementCriterion class). The subclasses implement the field.

    The measure of a cell is the largest absolute difference between the field value of the cell
    and that of any of its face neighbors, which are the cells that share a wall, or part of a
    wall, with it, using the field values averaged over the iterations of the averaging window of
    the refinement. A wall on the edge of the domain has no neighbors. The criterion flags a cell
    if its measure exceeds the \em maxChange threshold, with an excess equal to the ratio of the
    measure to the threshold, unless the field value of the cell is below the \em minValue floor.

    The measure, the threshold, and the floor are expressed in the internal units of SKIRT for the
    quantity of the field, because the units of a field selected by name are not known when the
    configuration is entered. Fields with different peak values therefore need different values for
    the threshold and the floor. By default, there is no floor. */
class GradientCriterion : public RefinementCriterion
{
    ITEM_ABSTRACT(GradientCriterion, RefinementCriterion, "a dynamic grid refinement criterion based on a gradient")

        PROPERTY_DOUBLE(maxChange, "the largest difference of the field value between neighboring cells")
        ATTRIBUTE_MIN_VALUE(maxChange, "]0")
        ATTRIBUTE_MAX_VALUE(maxChange, "1e300")

        PROPERTY_DOUBLE(minValue, "the field value below which a cell is never flagged")
        ATTRIBUTE_MIN_VALUE(minValue, "-1e300")
        ATTRIBUTE_MAX_VALUE(minValue, "1e300")
        ATTRIBUTE_DEFAULT_VALUE(minValue, "-1e300")
        ATTRIBUTE_DISPLAYED_IF(minValue, "Level2")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns true, because the test needs the face neighbors of a cell. */
    bool needsNeighbors() const override;

    /** This function returns the excess with which the criterion flags the spatial cell with index
        \f$m\f$, or zero if it does not flag the cell, as described in the class header. */
    double excess(int m, const CellWindow& window, const vector<int>& neighbors) const override;
};

////////////////////////////////////////////////////////////////////

#endif
