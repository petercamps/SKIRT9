/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VALUECRITERION_HPP
#define VALUECRITERION_HPP

#include "RefinementCriterion.hpp"

////////////////////////////////////////////////////////////////////

/** ValueCriterion is the abstract base class for dynamic grid refinement criteria that flag each
    spatial cell in which the field lies in a given range (see the RefinementCriterion class). The
    criterion does not look at the neighboring cells. The subclasses implement the field.

    A cell is flagged if its field value, averaged over the iterations of the averaging window of
    the refinement, lies in the range from \em minValue to \em maxValue, with an excess of one.
    Both bounds are inclusive, and are expressed in the internal units of SKIRT for the quantity of
    the field. Either bound can be left at its default value, which leaves the range open on that
    side.

    A value criterion expresses that a given kind of material must be resolved down to a given
    scale, wherever it occurs. Because child cells inherit their parent's state, a cell flagged by a
    value criterion remains flagged after it has been subdivided, at least until later iterations
    change its state. The final resolution is therefore determined by the limits of the criterion,
    and not by its test. Setup reports a fatal error if neither a level limit nor a column density
    limit applies. */
class ValueCriterion : public RefinementCriterion
{
    ITEM_ABSTRACT(ValueCriterion, RefinementCriterion, "a dynamic grid refinement criterion based on a value range")

        PROPERTY_DOUBLE(minValue, "the smallest field value for which a cell is flagged")
        ATTRIBUTE_MIN_VALUE(minValue, "-1e300")
        ATTRIBUTE_MAX_VALUE(minValue, "1e300")
        ATTRIBUTE_DEFAULT_VALUE(minValue, "-1e300")

        PROPERTY_DOUBLE(maxValue, "the largest field value for which a cell is flagged")
        ATTRIBUTE_MIN_VALUE(maxValue, "-1e300")
        ATTRIBUTE_MAX_VALUE(maxValue, "1e300")
        ATTRIBUTE_DEFAULT_VALUE(maxValue, "1e300")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the range is not empty, and that a level limit or a column
        density limit applies. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns false, because the test does not need the face neighbors of a cell.
        */
    bool needsNeighbors() const override;

    /** This function returns one if the field value of the spatial cell with index \f$m\f$ lies in
        the range of the criterion, and zero otherwise. */
    double excess(int m, const CellWindow& window, const vector<int>& neighbors) const override;
};

////////////////////////////////////////////////////////////////////

#endif
