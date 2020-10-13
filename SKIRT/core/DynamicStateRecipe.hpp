/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DYNAMICSTATERECIPE_HPP
#define DYNAMICSTATERECIPE_HPP

#include "Array.hpp"
#include "SimulationItem.hpp"
class MaterialState;

////////////////////////////////////////////////////////////////////

/** DynamicStateRecipe is an abstract class representing recipes for dynamically adjusting the
    medium state by iterating over the radiation field. The purpose of the dynamic medium state
    mechanism is to allow self-consistent calculation of medium state variables. A prototypical
    example is the radiative destruction of dust grains.

    For each iteration step, the simulation launches a set of photon packets to track the radiation
    field based on the current values of the medium state variables. At the end of the iteration
    step, the dynamic medium state recipes configured for the simulation are given a chance to
    update the medium state variables based on the established radiation field. The iteration
    continues until all recipes indicate convergence (or the configured maximum number of
    iterations has been reached).

    At the end of each iteration step, the functions for each configured recipe are called in the
    following order: the beginUpdate() function; the update() function for each spatial cell and
    for each medium component in the simulation; and finally the endUpdate() function. The latter
    function also returns a Boolean that indicates whether the iteration can be considered to have
    converged for this recipe.

    A recipe is allowed to discover information about the configuration of the simulation, for
    example to determine the type of material mixes associated with the various medium components.
    However, it should not update any information outside its own data members except through the
    material states passed to the update() function. Furthermore, the operation of the update()
    function should solely depend on and affect changes to the material state passed to that
    particular invocation of the function. In other words, each spatial cell/component must be
    treated independently from other spatial cells/components.

    The functions for multiple recipes are invoked in the user-configured order for each particular
    cell (but may be otherwise interleaved in any order). Assuming two recipes A and B, configured
    in this order, recipe B will thus always see the state as it has been adjusted by recipe A.

    These requirements and guarantees allow the update() function to be called for different
    spatial cell/component combinations in parallel execution threads and in multiple parallel
    processes, and thus in arbitrary order. A possible execution sequence could be:

    \verbatim
    for r in recipes:
        r.beginUpdate()

    parallel for m in spatial cells:
        for r in recipes:
            for h in medium components:
                r.update(m, h)

    for r in recipes:
        r.endUpdate()
    \endverbatim
    */
class DynamicStateRecipe : public SimulationItem
{
    ITEM_ABSTRACT(DynamicStateRecipe, SimulationItem, "a dynamic medium state recipe")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function is called before updating begins at the end of each iteration step. The \em
        numCells argument specifies the number of spatial cells in the simulation. The
        implementation in a subclass can perform the following tasks, if applicable:

        - gather and cache information about the simulation at large (this is easier done here than
        during setup because recipe instances are being setup before the medium system of which
        they are a part).

        - initialize the data members that track convergence data for the recipe in case the recipe
        requires more information than the number of updated spatial cells, which is automatically
        provided to the endUpdate() function. */
    virtual void beginUpdate(int numCells) = 0;

    /** This function is called during updating at the end of each iteration step. If needed, it
        updates the specified material state (corresponding to a particular spatial cell and medium
        component) based on the specified radiation field (corresponding to that same cell and
        discretized on the simulation's radiation field wavelength grid as returned by the
        Configuration::radiationFieldWLG() function). The function returns true if the material
        state has been updated and false if the material state has been left unchanged.

        If applicable, the implementation in a subclass should also update the data members that
        track additional convergence data for the recipe. Because the update() function can be
        called from multiple parallel execution threads, those updates must be performed atomically
        or be properly protected by a lock. */
    virtual bool update(MaterialState* state, const Array& Jv) = 0;

    /** This function is called after updating has completed at the end of each iteration step. It
        returns true if the recipe has converged and false if not. The \em numCells and \em
        numUpdated arguments specify respectively the number of spatial cells in the simulation and
        the number of cells updated during this update cycle.

        In case the recipe tracks additional information in the update() function, the tracked data
        must be aggregated across processes, cause the update() function can be called from
        separate parallel processes.

        In addition to determining the binary yes/no convergence result, the implementation in a
        subclass should also issue a log message reporting the degree of convergence. */
    virtual bool endUpdate(int numCells, int numUpdated) = 0;
};

////////////////////////////////////////////////////////////////////

#endif