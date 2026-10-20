/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TRIVIALSTAGEDGASMIX_HPP
#define TRIVIALSTAGEDGASMIX_HPP

#include "StateVariableList.hpp"
#include "TrivialGasMix.hpp"

////////////////////////////////////////////////////////////////////

/** The TrivialStagedGasMix class represents a trivial gas mix, with the optical properties of the
    TrivialGasMix class, that has a primary dynamic medium state updated in four stages, as
    described for the MaterialMix class. It is provided for testing the staged update protocol and
    serves no physical purpose: the dynamic state does not affect the optical properties.

    The update mimics a material mix that obtains its state from an external calculation, which is
    performed for the nodes of a parameter grid rather than for each spatial cell:

    - Stage 1 determines for each spatial cell with material the node on a grid in the logarithm
      of the number density, with the configured node spacing in dex, and stores the node index in
      the medium state of the cell. The node with index \f$k\f$ represents the number density
      \f$10^{k\Delta}\,\mathrm{m}^{-3}\f$, where \f$\Delta\f$ is the node spacing, and a cell with
      number density \f$n\f$ belongs to the node with index \f$\lfloor\log_{10}(n/\mathrm{m}^{-3})
      / \Delta\rfloor\f$.

    - Stage 2 collects the distinct node indices of all cells, and calculates the number density
      represented by each of these nodes, distributing the nodes over the available threads and
      processes. It keeps the resulting table in the material mix.

    - Stage 3 looks up the table entry for the node of each cell, and stores the number density of
      the node in the medium state of the cell, together with the number of nodes in the table.

    - Stage 4 logs the number of nodes in the table, and releases the table.

    The update does not use the radiation field, so that the dynamic state does not depend on the
    random numbers of the simulation, nor on the number of threads and processes. The number of
    nodes stored in each cell equals the number of distinct nodes in the complete spatial grid,
    which verifies that stage 2 sees the node indices determined by all processes in stage 1.

    The specific state variables are, in addition to the number density, the node index (named \c
    node), the number density of the node (named \c n_node), and the number of nodes in the table
    (named \c numNodes). */
class TrivialStagedGasMix : public TrivialGasMix
{
    ITEM_CONCRETE(TrivialStagedGasMix, TrivialGasMix,
                  "A trivial gas mix with a staged dynamic medium state update for testing purposes")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TrivialStagedGasMix, "Level3")
        ATTRIBUTE_TYPE_INSERT(TrivialStagedGasMix, "CustomMediumState,DynamicState")

        PROPERTY_DOUBLE(nodeSpacing, "the spacing of the grid nodes in the logarithm of the number density (dex)")
        ATTRIBUTE_MIN_VALUE(nodeSpacing, "]0")
        ATTRIBUTE_MAX_VALUE(nodeSpacing, "10]")
        ATTRIBUTE_DEFAULT_VALUE(nodeSpacing, "0.5")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function assembles the list of specific state variables. */
    void setupSelfBefore() override;

    //======== Capabilities =======

public:
    /** This function returns DynamicStateType::Primary, indicating that the medium state of this
        material mix is updated after each primary emission segment. */
    DynamicStateType hasDynamicMediumState() const override;

    /** This function returns true, indicating that the medium state of this material mix is
        updated in stages. */
    bool hasStagedStateUpdate() const override;

    //======== Medium state setup =======

public:
    /** This function returns the list of specific state variables used by this material mix: the
        number density, the node index, the number density of the node, and the number of nodes in
        the table. */
    const StateVariableList& specificStateVariableInfo() const override;

    //======== Medium state updates =======

public:
    /** This function implements stage 1 of the update: for a spatial cell with material, it stores
        the index of the node for the number density of the cell, and returns an updated status. */
    UpdateStatus prepareSpecificStateUpdate(MaterialState* state, const Array& Jv) const override;

    /** This function implements stage 2 of the update: it collects the distinct node indices of
        the spatial cells with material for the specified medium component, and builds the table
        with the number density of each node, distributing the calculation over threads and
        processes. */
    void beginSpecificStateUpdate(const MediumState& state, int h) const override;

    /** This function implements stage 3 of the update: for a spatial cell with material, it stores
        the number density of the node of the cell, found in the table built by stage 2, and the
        number of nodes in the table, and returns an updated status. The function throws a fatal
        error if the node is missing from the table. */
    UpdateStatus updateSpecificState(MaterialState* state, const Array& Jv) const override;

    /** This function implements stage 4 of the update: it logs the number of nodes in the table
        built by stage 2, and releases the table. */
    void endSpecificStateUpdate(const MediumState& state, int h) const override;

    //======================== Data Members ========================

private:
    // the specific state variables and their handles - initialized during setup
    StateVariableList _variables;
    StateScalar _node;
    StateScalar _nodeDensity;
    StateScalar _numNodes;

    // the table built by stage 2 and released by stage 4: the sorted node indices and the corresponding number densities
    mutable vector<int> _tableNodes;
    mutable Array _tableDensities;
};

////////////////////////////////////////////////////////////////////

#endif
