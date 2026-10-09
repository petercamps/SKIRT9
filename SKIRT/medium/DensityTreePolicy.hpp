/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DENSITYTREEPOLICY_HPP
#define DENSITYTREEPOLICY_HPP

#include "MaterialMix.hpp"
#include "TreePolicy.hpp"

//////////////////////////////////////////////////////////////////////

/** DensityTreePolicy is a tree subdivision policy that limits the fraction of the material of a
    given type (dust, electrons, or gas) contained in each cell. A node is subdivided if it contains
    a fraction of the total amount of material of the selected type that exceeds the configured
    maximum fraction \f$\delta_\text{max}\f$.

    For dust, the policy uses mass and mass density, which is the appropriate quantity in case
    multiple dust medium components have a different mass per hydrogen atom. For electrons and gas,
    it uses number and number density. The description below refers to mass and mass density; for
    electrons and gas, read number and number density instead.

    The total mass in the model, \f$M_\text{model}\f$, and the mass \f$M\f$ inside a node are
    obtained by summing the corresponding quantity for each medium component of the selected
    material type. The mass of a medium component inside the node is calculated exactly if the
    component offers the MassInBoxInterface, and is otherwise estimated from density samples (see
    the TreeNodeEvaluation class). The fraction of the mass within the node is then \f[\delta =
    \frac{M}{M_\text{model}}.\f]

    If the simulation has no medium components of the selected material type, setup reports a fatal
    error. */
class DensityTreePolicy : public TreePolicy
{
    /** The enumeration type indicating the material type to which the criterion applies. */
    ENUM_DEF(MaterialType, Dust, Electrons, Gas)
        ENUM_VAL(MaterialType, Dust, "dust")
        ENUM_VAL(MaterialType, Electrons, "electrons")
        ENUM_VAL(MaterialType, Gas, "gas")
    ENUM_END()

    ITEM_CONCRETE(DensityTreePolicy, TreePolicy,
                  "a tree subdivision policy limiting the fraction of material in a cell")

        PROPERTY_ENUM(materialType, MaterialType, "the material type to which the criterion applies")
        ATTRIBUTE_DEFAULT_VALUE(materialType, "DustMix:Dust;GasMix:Gas;Electrons")

        PROPERTY_DOUBLE(maxFraction, "the maximum fraction of the material contained in each cell")
        ATTRIBUTE_MIN_VALUE(maxFraction, "]0")
        ATTRIBUTE_MAX_VALUE(maxFraction, "1e-2]")
        ATTRIBUTE_DEFAULT_VALUE(maxFraction, "1e-6")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function obtains and caches information used by the needsSubdivide() function. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns true if the specified node contains a fraction of the material of
        the selected type that exceeds the configured maximum fraction, and false otherwise. */
    bool needsSubdivide(TreeNodeEvaluation& node) const override;

    //======================== Data Members ========================

private:
    // data members initialized by setupSelfBefore()
    MaterialMix::MaterialType _type{MaterialMix::MaterialType::Dust};  // the selected material type
    double _total{0.};  // the total amount of material of that type in the model
};

//////////////////////////////////////////////////////////////////////

#endif
