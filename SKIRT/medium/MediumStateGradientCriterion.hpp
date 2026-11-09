/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEDIUMSTATEGRADIENTCRITERION_HPP
#define MEDIUMSTATEGRADIENTCRITERION_HPP

#include "GradientCriterion.hpp"
#include "StateVariableSelection.hpp"
class DynamicStateRecipe;
class MediumSystem;

////////////////////////////////////////////////////////////////////

/** A MediumStateGradientCriterion flags the spatial cells in which a medium state variable differs
    strongly from that in a neighboring cell, as described for the GradientCriterion class. It is
    used for dynamic grid refinement (see the DynamicRefinementOptions class).

    The variable is selected by its short name, as declared by the material mix that offers it, for
    example \c x_HI for the neutral hydrogen fraction of a DiffuseIonizedGasMix, which resolves the
    ionization front. It is taken from the first medium component whose material mix offers a
    scalar custom state variable with that name. Setup reports a fatal error listing the available
    names if no component offers the requested one, and also if that component has no dynamic
    medium state, because the variable then does not change between iterations. The criterion takes
    part in the iteration loops in which the variable changes: the primary and merged loops for a
    component with a primary dynamic medium state, and the secondary and merged loops for a
    component with a secondary dynamic medium state. If the simulation has dynamic medium state
    recipes, which may update the state of any component, the criterion also takes part in the
    primary and merged loops.

    The noise test of the \em significance property (see the GradientCriterion class) needs an
    estimate of the noise of the variable. The criterion obtains it from the first dynamic medium
    state recipe that sets the variable and can recalculate it for a scaled radiation field (see
    the DynamicStateRecipe::offersScaledFieldStateVariable() function), such as a
    DustDestructionRecipe for the dynamic density fraction \c fdyn_f of a fragmented dust mix. The
    noise of the variable in a cell is then the difference between its values for the radiation
    field of the cell multiplied by \f$1+k\sigma\f$ and by \f$1-k\sigma\f$, divided by \f$2k\f$,
    where \f$k\f$ is the significance and \f$\sigma\f$ the relative noise of the radiation field
    in the cell (see the MediumSystem::radiationFieldNoise() function). For a variable that varies
    smoothly with the strength of the radiation field, this is the standard error of the variable;
    for a variable that jumps, such as the density fraction at a sharp destruction front, it
    accounts for a jump within \f$k\f$ standard errors of the radiation field, so that the
    refinement of the front stops where its cells are too noisy to resolve it. If no recipe offers
    the variable, the criterion ignores the \em significance property: setup logs a warning if the
    property is nonzero, and the medium system does not track the noise for this criterion.

    The number density for the column density limit is the one stored in the medium state for the
    medium component that provides the variable. Its meaning depends on the material mix. For the
    DiffuseIonizedGasMix, it counts hydrogen plus helium nuclei, so that a limit intended as a
    hydrogen column density must be multiplied by one plus the helium abundance. */
class MediumStateGradientCriterion : public GradientCriterion
{
    ITEM_CONCRETE(MediumStateGradientCriterion, GradientCriterion,
                  "a criterion on the difference of a state variable between neighboring cells")
        ATTRIBUTE_TYPE_ALLOWED_IF(MediumStateGradientCriterion, "CustomMediumState")

        PROPERTY_STRING(variable, "the name of the medium state variable")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function selects the medium state variable and, if the \em significance property is
        nonzero, the dynamic medium state recipe that provides the noise of the variable, logging a
        warning if there is none. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns true if the variable changes in the specified iteration loop. */
    bool changesIn(IterationHistory::Loop loop) const override;

    /** This function returns the current value of the variable in the spatial cell with index
        \f$m\f$. */
    double value(int m) const override;

    /** This function returns the number density of the medium component that provides the
        variable in the spatial cell with index \f$m\f$. */
    double numberDensity(int m) const override;

    /** This function returns true if the \em significance property is nonzero and a dynamic medium
        state recipe provides the noise of the variable, and false otherwise (see the class
        header). */
    bool usesNoise() const override;

    /** This function returns true if the criterion uses the noise of the variable, because the
        recipe then recalculates the variable twice for each cell, and false otherwise. */
    bool hasCostlyValue() const override;

    /** This function returns an estimate of the standard error of the variable in the spatial cell
        with index \f$m\f$, as described in the class header. A cell without any contributions to
        one of the noise batches of the radiation field has infinite noise. */
    double noise(int m) const override;

    /** This function returns a description of the variable. */
    string fieldDescription() const override;

    //======================== Data Members ========================

private:
    StateVariableSelection _selection;
    const MediumSystem* _ms{nullptr};
    const DynamicStateRecipe* _recipe{nullptr};  // the recipe that provides the noise, or nullptr if none
};

////////////////////////////////////////////////////////////////////

#endif
