/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEDIUMSTATECONDITION_HPP
#define MEDIUMSTATECONDITION_HPP

#include "RefinementCondition.hpp"
#include "StateVariableSelection.hpp"

////////////////////////////////////////////////////////////////////

/** A MediumStateCondition restricts a dynamic grid refinement criterion to the spatial cells in
    which a medium state variable lies in a given range (see the RefinementCriterion class). The
    variable is selected by its short name, as declared by the material mix that offers it, for
    example \c x_HI for the neutral hydrogen fraction of a DiffuseIonizedGasMix. It is taken from
    the first medium component whose material mix offers a scalar custom state variable with that
    name. Setup reports a fatal error listing the available names if no component offers the
    requested one. The variable of a condition may be static, i.e. it does not need to change
    between iterations.

    The condition holds in a cell if the value of the variable, averaged over the iterations of the
    averaging window of the refinement, lies in the range from \em minValue to \em maxValue. Both
    bounds are inclusive, and are expressed in the internal units of SKIRT for the quantity of the
    variable, because the units of a variable selected by name are not known when the
    configuration is entered. Either bound can be left at its default value, which leaves the range
    open on that side. */
class MediumStateCondition : public RefinementCondition
{
    ITEM_CONCRETE(MediumStateCondition, RefinementCondition,
                  "a condition requiring a medium state variable to lie in a given range")
        ATTRIBUTE_TYPE_ALLOWED_IF(MediumStateCondition, "CustomMediumState")

        PROPERTY_STRING(variable, "the name of the medium state variable")

        PROPERTY_DOUBLE(minValue, "the smallest value of the variable for which the condition holds")
        ATTRIBUTE_MIN_VALUE(minValue, "-1e300")
        ATTRIBUTE_MAX_VALUE(minValue, "1e300")
        ATTRIBUTE_DEFAULT_VALUE(minValue, "-1e300")

        PROPERTY_DOUBLE(maxValue, "the largest value of the variable for which the condition holds")
        ATTRIBUTE_MIN_VALUE(maxValue, "-1e300")
        ATTRIBUTE_MAX_VALUE(maxValue, "1e300")
        ATTRIBUTE_DEFAULT_VALUE(maxValue, "1e300")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function selects the medium state variable, and verifies that the range is not
        empty. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the current value of the variable in the spatial cell with index
        \f$m\f$. */
    double value(int m) const override;

    /** This function returns true if the specified value lies in the range of the condition. */
    bool holds(double value) const override;

    /** This function returns a description of the variable. */
    string description() const override;

    //======================== Data Members ========================

private:
    StateVariableSelection _selection;
};

////////////////////////////////////////////////////////////////////

#endif
