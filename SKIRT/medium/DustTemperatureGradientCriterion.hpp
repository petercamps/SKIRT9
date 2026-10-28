/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTTEMPERATUREGRADIENTCRITERION_HPP
#define DUSTTEMPERATUREGRADIENTCRITERION_HPP

#include "GradientCriterion.hpp"
class MediumSystem;

////////////////////////////////////////////////////////////////////

/** A DustTemperatureGradientCriterion flags the spatial cells in which the indicative dust
    temperature differs strongly from that in a neighboring cell, as described for the
    GradientCriterion class. It is used for dynamic grid refinement (see the
    DynamicRefinementOptions class), for example to resolve the structure near embedded sources or
    at the edge of an optically thick region.

    The indicative dust temperature of a cell is calculated by the
    MediumSystem::indicativeDustTemperature() function, which solves the energy balance equation
    for a representative grain of each dust component in the local radiation field, and averages
    the result over the dust components weighted by their mass in the cell. It is not stored in the
    medium state, but calculated from the radiation field when needed. The threshold and the floor
    of the criterion are expressed in kelvin. Setup reports a fatal error if the simulation has no
    dust or does not store the radiation field.

    The criterion takes part in all iteration loops, because each loop recalculates at least part of
    the radiation field. In the secondary emission loop, however, the primary radiation field is not
    recalculated, and a new cell receives a uniform share of its parent's primary radiation field.
    For a temperature that depends mostly on the primary radiation field, the jump in value at the
    boundary of a subdivided cell therefore does not diminish with further subdivision, so that the
    cells along that boundary may keep asking for subdivision until a limit or the cell cap stops
    them.

    The number density for the column density limit is the sum of the number densities of the dust
    medium components. */
class DustTemperatureGradientCriterion : public GradientCriterion
{
    ITEM_CONCRETE(DustTemperatureGradientCriterion, GradientCriterion,
                  "a criterion on the difference of the indicative dust temperature between neighboring cells")
        ATTRIBUTE_TYPE_ALLOWED_IF(DustTemperatureGradientCriterion, "DustMix&RadiationField")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the simulation has dust and stores the radiation field. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns true for any iteration loop. */
    bool changesIn(IterationHistory::Loop loop) const override;

    /** This function returns the indicative dust temperature in the spatial cell with index
        \f$m\f$. */
    double value(int m) const override;

    /** This function returns true, because calculating the indicative dust temperature is
        expensive. */
    bool hasCostlyValue() const override;

    /** This function returns the sum of the number densities of the dust medium components in the
        spatial cell with index \f$m\f$. */
    double numberDensity(int m) const override;

    /** This function returns a description of the field. */
    string fieldDescription() const override;

    //======================== Data Members ========================

private:
    const MediumSystem* _ms{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
