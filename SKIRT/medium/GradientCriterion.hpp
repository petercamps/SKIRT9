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
    the threshold and the floor. By default, there is no floor.

    If the \em logarithmic flag is set, the measure is instead the largest absolute difference
    between the natural logarithms of the field values of the cell and of its face neighbors. For
    small differences, this is the relative difference, so that the threshold is dimensionless (a
    value of 0.01 allows a change of about one percent between neighboring cells) and serves fields
    with any peak value, such as a temperature that varies by a large factor across the domain or a
    mean intensity that varies by orders of magnitude. The floor is still expressed in the units of
    the field. A cell with a field value that is not positive is then never flagged, and a neighbor
    with such a value is ignored.

    <b>Noise</b>

    The field values carry Monte Carlo noise, which grows as cells get smaller, because fewer
    photon packets cross a smaller cell. Noise can push the difference between neighboring cells
    over the threshold, in particular for the many cells whose difference is just below it; the
    children of a cell subdivided for that reason have smaller differences but are noisier still,
    so that the refinement may cascade to ever smaller cells. If the \em significance property is
    nonzero, the difference with each face neighbor is therefore reduced by \em significance times
    its estimated standard error, which combines the noise of the window-averaged field values of
    both cells, before it is compared with the threshold. A cell is thus flagged only if its
    difference exceeds the threshold with a confidence set by the significance factor. The noise of
    a field value is estimated by the subclass from the noise of the radiation field in the cell
    (see the MediumSystem::radiationFieldNoise() function); a subclass that provides no estimate is
    not affected. Where the noise is small compared to the threshold, the test changes little;
    where it is not, the effective threshold rises, so that a cell that is too noisy to resolve a
    gradient is not subdivided further. */
class GradientCriterion : public RefinementCriterion
{
    ITEM_ABSTRACT(GradientCriterion, RefinementCriterion, "a dynamic grid refinement criterion based on a gradient")

        PROPERTY_BOOL(logarithmic, "compare the logarithms of the field values rather than the values themselves")
        ATTRIBUTE_DEFAULT_VALUE(logarithmic, "false")

        PROPERTY_DOUBLE(maxChange, "the largest difference of the field value between neighboring cells")
        ATTRIBUTE_MIN_VALUE(maxChange, "]0")
        ATTRIBUTE_MAX_VALUE(maxChange, "1e300")

        PROPERTY_DOUBLE(significance, "the number of standard errors by which a difference must exceed the "
                                      "threshold (zero means no noise test)")
        ATTRIBUTE_MIN_VALUE(significance, "[0")
        ATTRIBUTE_MAX_VALUE(significance, "100]")
        ATTRIBUTE_DEFAULT_VALUE(significance, "0")
        ATTRIBUTE_DISPLAYED_IF(significance, "Level2")

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

    /** This function returns true if the \em significance property is nonzero, and false
        otherwise. */
    bool usesNoise() const override;

    /** This function returns the excess with which the criterion flags the spatial cell with index
        \f$m\f$, or zero if it does not flag the cell, as described in the class header. */
    double excess(int m, const CellWindow& window, const CellWindow* noiseWindow,
                  const vector<int>& neighbors) const override;
};

////////////////////////////////////////////////////////////////////

#endif
