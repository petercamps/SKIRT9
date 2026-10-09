/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DISPERSIONTREEPOLICY_HPP
#define DISPERSIONTREEPOLICY_HPP

#include "MaterialMix.hpp"
#include "TreePolicy.hpp"

//////////////////////////////////////////////////////////////////////

/** DispersionTreePolicy is a tree subdivision policy that limits the dispersion of the density of
    a given material type (dust, electrons, or gas) within each cell. A node is subdivided if the
    density dispersion measure \f$q\f$ calculated for the node exceeds the configured maximum
    dispersion \f$q_\text{max}\f$.

    The density at a given position is obtained by summing the density of each medium component of
    the selected material type; it is the mass density for dust and the number density for electrons
    and gas. The density is sampled in \f$N\f$ random positions distributed uniformly across the
    volume of the node, where \f$N\f$ is configured in the SamplingOptions of the medium system (see
    the TreeNodeEvaluation class). The dispersion measure \f$q\f$ is then determined as \f[ q =
    \begin{cases} \;\dfrac{\rho_{\text{max}}-\rho_{\text{min}}}{\rho_{\text{max}}} & \quad\text{if
    $\rho_{\text{max}}>0$,} \\ \;0 & \quad\text{if $\rho_{\text{max}}=0$,} \end{cases} \f] where
    \f$\rho_{\text{min}}\f$ and \f$\rho_{\text{max}}\f$ are the smallest and largest sampled density
    values. The quantity \f$q\f$ is a simple measure for the uniformity of the density within the
    node: for a constant density, \f$q=0\f$, whereas \f$q\f$ approaches 1 if a steep gradient is
    present. The special case \f$\rho_{\text{max}}=0\f$ covers an empty node, for which the uniform
    value \f$q=0\f$ is returned.

    Nodes that contain a sharp edge with empty space on one side have \f$q=1\f$, so that they will
    continue to be subdivided for ever if \f$q_\text{max}<1\f$. It is thus important to always
    configure a reasonable maximum subdivision level for the grid when using this policy.

    If the simulation has no medium components of the selected material type, setup reports a fatal
    error. */
class DispersionTreePolicy : public TreePolicy
{
    /** The enumeration type indicating the material type to which the criterion applies. */
    ENUM_DEF(MaterialType, Dust, Electrons, Gas)
        ENUM_VAL(MaterialType, Dust, "dust")
        ENUM_VAL(MaterialType, Electrons, "electrons")
        ENUM_VAL(MaterialType, Gas, "gas")
    ENUM_END()

    ITEM_CONCRETE(DispersionTreePolicy, TreePolicy,
                  "a tree subdivision policy limiting the density dispersion in a cell")
        ATTRIBUTE_TYPE_DISPLAYED_IF(DispersionTreePolicy, "Level2")

        PROPERTY_ENUM(materialType, MaterialType, "the material type to which the criterion applies")
        ATTRIBUTE_DEFAULT_VALUE(materialType, "DustMix:Dust;GasMix:Gas;Electrons")

        PROPERTY_DOUBLE(maxDispersion, "the maximum density dispersion in each cell")
        ATTRIBUTE_MIN_VALUE(maxDispersion, "]0")
        ATTRIBUTE_MAX_VALUE(maxDispersion, "1]")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function obtains and caches information used by the needsSubdivide() function. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns true if the density dispersion measure for the specified node
        exceeds the configured maximum dispersion, and false otherwise. */
    bool needsSubdivide(TreeNodeEvaluation& node) const override;

    //======================== Data Members ========================

private:
    // data members initialized by setupSelfBefore()
    MaterialMix::MaterialType _type{MaterialMix::MaterialType::Dust};  // the selected material type
};

//////////////////////////////////////////////////////////////////////

#endif
