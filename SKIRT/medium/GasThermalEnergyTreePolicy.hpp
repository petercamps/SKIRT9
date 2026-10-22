/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GASTHERMALENERGYTREEPOLICY_HPP
#define GASTHERMALENERGYTREEPOLICY_HPP

#include "TreePolicy.hpp"
class Medium;

//////////////////////////////////////////////////////////////////////

/** GasThermalEnergyTreePolicy is a tree subdivision policy that limits the fraction of the total
    thermal energy of the gas contained in each cell. It mirrors the DensityTreePolicy for gas, but
    weighs the number density of the gas with its temperature. A node is subdivided if its thermal
    energy, \f[ E = \frac{3}{2}\,k_\text{B} \int_V \sum_h n_h({\bf{r}})\, T_h({\bf{r}})\,
    {\text{d}}V, \f] exceeds the configured fraction of the total thermal energy of the gas in the
    spatial domain, where \f$n_h\f$ and \f$T_h\f$ are the number density and temperature of gas
    medium component \f$h\f$. The constant factor cancels in this fraction and is not actually
    calculated.

    The number density \f$n_h\f$ counts the entities of the material mix of each medium component,
    which differ between mixes: for example hydrogen atoms, neutral hydrogen atoms, or the
    molecules of a trace species. Neither these differences nor the number of free particles per
    entity (such as the electrons and helium atoms in ionized gas) are taken into account. The
    policy is therefore intended for a single gas medium component, or for gas medium components
    whose entities are comparable.

    The tree is constructed before the medium state, including the gas temperature, has been
    initialized. The policy therefore uses the temperature imported by the medium components, i.e.
    the temperature column of imported media whose material mix or material mix family requests the
    temperature as an imported quantity (see the Medium::parameterInfo() and Medium::parameters()
    functions). Gas medium components that do not import a temperature, including all geometric
    media, are ignored, and setup logs a warning listing them. If none of the gas medium components
    imports a temperature, setup reports a fatal error.

    The thermal energy of a medium component in a node is estimated as the number of entities in the
    node multiplied by the average temperature in the node. The number of entities is calculated
    exactly if the medium component offers the MassInBoxInterface, and is otherwise estimated from
    the number density sampled at the random positions shared by the policies evaluating the node
    (see the TreeNodeEvaluation class). A particle medium with a material mix family whose members
    may have a different mass per entity (see Medium::hasUniformMass()) does not offer this
    interface, so that the accuracy of the estimate then depends on the number of density samples
    configured for the medium system. The average temperature is the average of the
    temperature sampled at the same positions, weighted by the sampled number density. If none of
    the samples hits any material while the node does contain material (which can happen for a small
    particle in a large node), the average temperature of the medium component is used instead. This
    average temperature is estimated during setup as the average over a large number of positions
    drawn from the spatial distribution of the component, and the total thermal energy of the
    component as its total number of entities multiplied by its average temperature. Without the
    MassInBoxInterface, the estimate for a node reduces to \f$E \propto (V/N) \sum_i n({\bf{r}}_i)\,
    T({\bf{r}}_i)\f$, with \f$N\f$ samples at positions \f${\bf{r}}_i\f$ in a node with volume
    \f$V\f$. */
class GasThermalEnergyTreePolicy : public TreePolicy
{
    ITEM_CONCRETE(GasThermalEnergyTreePolicy, TreePolicy,
                  "a tree subdivision policy limiting the fraction of the gas thermal energy in a cell")
        ATTRIBUTE_TYPE_ALLOWED_IF(GasThermalEnergyTreePolicy, "GasMix")
        ATTRIBUTE_TYPE_DISPLAYED_IF(GasThermalEnergyTreePolicy, "Level2")

        PROPERTY_DOUBLE(maxFraction, "the maximum fraction of the gas thermal energy contained in each cell")
        ATTRIBUTE_MIN_VALUE(maxFraction, "]0")
        ATTRIBUTE_MAX_VALUE(maxFraction, "1e-2]")
        ATTRIBUTE_DEFAULT_VALUE(maxFraction, "1e-6")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function selects the gas medium components that import a temperature, and estimates
        their total thermal energy. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns true if the estimated thermal energy of the gas in the specified node
        exceeds the configured fraction of the total thermal energy, and false otherwise. */
    bool needsSubdivide(TreeNodeEvaluation& node) const override;

    //======================== Data Members ========================

private:
    // data members initialized by setupSelfBefore()
    vector<int> _hv;            // the indices of the gas medium components that import a temperature
    vector<Medium*> _media;     // the corresponding medium components
    vector<int> _indexv;        // the index of the temperature among the imported parameters of each component
    vector<double> _averageTv;  // the average temperature of each of these medium components
    double _maxEnergy{0.};      // the maximum thermal energy in a node, without the constant factor
};

//////////////////////////////////////////////////////////////////////

#endif
