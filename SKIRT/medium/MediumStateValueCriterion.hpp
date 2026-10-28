/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEDIUMSTATEVALUECRITERION_HPP
#define MEDIUMSTATEVALUECRITERION_HPP

#include "StateVariableSelection.hpp"
#include "ValueCriterion.hpp"

////////////////////////////////////////////////////////////////////

/** A MediumStateValueCriterion flags the spatial cells in which a medium state variable lies in a
    given range, as described for the ValueCriterion class. It is used for dynamic grid refinement
    (see the DynamicRefinementOptions class), for example to resolve the layer in which a given ion
    is concentrated down to a given column density per cell. Setup reports a fatal error if neither
    a level limit nor a column density limit applies to the criterion.

    The variable is selected by its short name, as declared by the material mix that offers it, for
    example \c x_NII for the fraction of singly ionized nitrogen of a DiffuseIonizedGasMix. It is
    taken from the first medium component whose material mix offers a scalar custom state variable
    with that name. Setup reports a fatal error listing the available names if no component offers
    the requested one, and also if that component has no dynamic medium state, because the variable
    then does not change between iterations. The criterion takes part in the iteration loops in
    which the variable changes: the primary and merged loops for a component with a primary dynamic
    medium state, and the secondary and merged loops for a component with a secondary dynamic medium
    state. If the simulation has dynamic medium state recipes, which may update the state of any
    component, the criterion also takes part in the primary and merged loops.

    The number density for the column density limit is the one stored in the medium state for the
    medium component that provides the variable. Its meaning depends on the material mix. For the
    DiffuseIonizedGasMix, it counts hydrogen plus helium nuclei, so that a limit intended as a
    hydrogen column density must be multiplied by one plus the helium abundance. */
class MediumStateValueCriterion : public ValueCriterion
{
    ITEM_CONCRETE(MediumStateValueCriterion, ValueCriterion, "a criterion on the value of a medium state variable")
        ATTRIBUTE_TYPE_ALLOWED_IF(MediumStateValueCriterion, "CustomMediumState")

        PROPERTY_STRING(variable, "the name of the medium state variable")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function selects the medium state variable. */
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

    /** This function returns a description of the variable. */
    string fieldDescription() const override;

    //======================== Data Members ========================

private:
    StateVariableSelection _selection;
};

////////////////////////////////////////////////////////////////////

#endif
