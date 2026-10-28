/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TRIVIALBOXCRITERION_HPP
#define TRIVIALBOXCRITERION_HPP

#include "Box.hpp"
#include "ValueCriterion.hpp"
class MediumSystem;
class TreeSpatialGrid;

////////////////////////////////////////////////////////////////////

/** The TrivialBoxCriterion class is a dynamic grid refinement criterion intended for testing
    purposes only. Its field is one for each spatial cell whose center lies inside the configured
    box, and zero for all other cells, so that, as a value criterion with a minimum value of one
    (or any value between zero and one), it flags the cells inside the box (see the ValueCriterion
    class). Because the field depends only on the geometry of the grid, the refinement decisions
    do not depend on the radiation field or on the random numbers of the simulation, which allows
    testing the mechanics of dynamic refinement, including the conditions, the limits, and the
    selection of the candidates under the cell cap, without Monte Carlo noise.

    The criterion takes part in all iteration loops. The number density for the column density
    limit is the sum of the number densities of all medium components. */
class TrivialBoxCriterion : public ValueCriterion
{
    ITEM_CONCRETE(TrivialBoxCriterion, ValueCriterion, "a criterion flagging the cells in a box, for testing purposes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TrivialBoxCriterion, "Level3")

        PROPERTY_DOUBLE(minX, "the start point of the box in the X direction")
        ATTRIBUTE_QUANTITY(minX, "length")

        PROPERTY_DOUBLE(maxX, "the end point of the box in the X direction")
        ATTRIBUTE_QUANTITY(maxX, "length")

        PROPERTY_DOUBLE(minY, "the start point of the box in the Y direction")
        ATTRIBUTE_QUANTITY(minY, "length")

        PROPERTY_DOUBLE(maxY, "the end point of the box in the Y direction")
        ATTRIBUTE_QUANTITY(maxY, "length")

        PROPERTY_DOUBLE(minZ, "the start point of the box in the Z direction")
        ATTRIBUTE_QUANTITY(minZ, "length")

        PROPERTY_DOUBLE(maxZ, "the end point of the box in the Z direction")
        ATTRIBUTE_QUANTITY(maxZ, "length")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the box is not empty. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns true for any iteration loop. */
    bool changesIn(IterationHistory::Loop loop) const override;

    /** This function returns one if the center of the spatial cell with index \f$m\f$ lies inside
        the box, and zero otherwise. */
    double value(int m) const override;

    /** This function returns the sum of the number densities of all medium components in the
        spatial cell with index \f$m\f$. */
    double numberDensity(int m) const override;

    /** This function returns a description of the field. */
    string fieldDescription() const override;

    //======================== Data Members ========================

private:
    const MediumSystem* _ms{nullptr};
    const TreeSpatialGrid* _grid{nullptr};
    Box _box;
};

////////////////////////////////////////////////////////////////////

#endif
