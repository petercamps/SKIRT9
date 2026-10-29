/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTLUMINOSITYCRITERION_HPP
#define DUSTLUMINOSITYCRITERION_HPP

#include "RefinementCriterion.hpp"
class MediumSystem;

////////////////////////////////////////////////////////////////////

/** A DustLuminosityCriterion flags the spatial cells that absorb more than a given fraction of the
    total luminosity absorbed by dust, so that no cell emits more than that fraction of the dust
    emission. It is used for dynamic grid refinement (see the DynamicRefinementOptions class) in
    simulations with dust emission, to resolve the regions that dominate the dust emission, such as
    the dust near embedded sources or the self-absorbing core of an optically thick galaxy, without
    resolving dust that hardly emits. It is the counterpart for the absorbed luminosity of the
    DensityTreePolicy, which limits the fraction of the dust mass in each cell when the grid is
    constructed.

    The field of the criterion is the fraction of the total dust-absorbed luminosity that is
    absorbed in the cell, which includes the luminosity absorbed during primary and secondary
    emission, as accumulated at the exact wavelength of each photon packet (see the MediumSystem
    class). If the \em primaryOnly flag is set, the field is instead the fraction of the
    dust-absorbed primary luminosity, so that the criterion resolves the regions where the radiation
    of the primary sources is absorbed, without subdividing the cells that absorb mostly the
    radiation emitted by the dust itself. A cell is flagged if this fraction, averaged over the
    iterations of the averaging window, exceeds the \em maxFraction threshold, with an excess equal
    to the ratio of the fraction to the threshold, so that the cells with the largest luminosity are
    subdivided first when the cell cap applies.

    Because a subdivided cell divides its absorbed luminosity among its children, the criterion
    limits itself: it stops flagging a cell once the cell absorbs less than the threshold, and it
    therefore needs no level limit or column density limit. In an optically thick cell, however,
    most of the luminosity is absorbed by the children facing the radiation, so that these continue
    to be subdivided, and around a bright source embedded in very opaque material the refinement
    follows the absorbing layer to ever smaller cells. In addition, Monte Carlo noise lets cells
    close to the threshold cross it in later averaging windows. A level limit or the \em maxRounds
    option of the refinement therefore bounds the refinement in practice.

    Without the \em primaryOnly flag, the criterion takes part in all iteration loops, because each
    loop recalculates at least part of the absorbed luminosity. In the primary emission loop, only
    the primary luminosity is absorbed, and in the secondary emission loop, the primary luminosity
    of a new cell is a share of that of its parent in proportion to its volume, because the primary
    radiation field is not recalculated. With the \em primaryOnly flag, the criterion takes part
    only in the primary and merged loops, which recalculate the primary radiation field; a
    simulation that iterates over primary emission before iterating over secondary emission thus
    refines the grid during the primary iterations, and keeps it during the secondary iterations.
    The number density for the column density limit is the sum of the number densities of the dust
    medium components. Setup reports a fatal error if the simulation does not accumulate the
    luminosity absorbed by dust, which requires dust emission and a dust medium. */
class DustLuminosityCriterion : public RefinementCriterion
{
    ITEM_CONCRETE(DustLuminosityCriterion, RefinementCriterion,
                  "a criterion on the fraction of the dust-absorbed luminosity in a cell")
        ATTRIBUTE_TYPE_ALLOWED_IF(DustLuminosityCriterion, "DustMix&DustEmission")

        PROPERTY_DOUBLE(maxFraction, "the maximum fraction of the dust-absorbed luminosity in each cell")
        ATTRIBUTE_MIN_VALUE(maxFraction, "]0")
        ATTRIBUTE_MAX_VALUE(maxFraction, "1]")
        ATTRIBUTE_DEFAULT_VALUE(maxFraction, "1e-5")

        PROPERTY_BOOL(primaryOnly, "consider only the luminosity absorbed from primary sources")
        ATTRIBUTE_DEFAULT_VALUE(primaryOnly, "false")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the simulation has dust emission and a dust medium. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns true for the loops in which the field changes, as described in the
        class header. */
    bool changesIn(IterationHistory::Loop loop) const override;

    /** This function obtains the total luminosity absorbed by dust, which normalizes the field. */
    void prepareValues() const override;

    /** This function returns the fraction of the total dust-absorbed luminosity that is absorbed in
        the spatial cell with index \f$m\f$. */
    double value(int m) const override;

    /** This function returns the sum of the number densities of the dust medium components in the
        spatial cell with index \f$m\f$. */
    double numberDensity(int m) const override;

    /** This function returns a description of the field. */
    string fieldDescription() const override;

    /** This function returns false, because the test does not need the face neighbors of a cell. */
    bool needsNeighbors() const override;

    /** This function returns the ratio of the window-averaged field value of the spatial cell with
        index \f$m\f$ to the \em maxFraction threshold if the ratio exceeds one, and zero
        otherwise. */
    double excess(int m, const CellWindow& window, const vector<int>& neighbors) const override;

    //======================== Data Members ========================

private:
    const MediumSystem* _ms{nullptr};
    mutable double _total{0.};  // the total dust-absorbed luminosity, obtained by prepareValues()
};

////////////////////////////////////////////////////////////////////

#endif
