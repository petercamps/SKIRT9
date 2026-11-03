/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTDESTRUCTIONRECIPE_HPP
#define DUSTDESTRUCTIONRECIPE_HPP

#include "DynamicStateRecipe.hpp"
class FragmentDustMixDecorator;
class MediumSystem;

////////////////////////////////////////////////////////////////////

/** DustDestructionRecipe is an abstract class representing recipes for dynamically calculating
    radiative destruction of dust grains by iterating over the radiation field. The class derives
    from DynamicStateRecipe so that it fits in the overall dynamic medium state framework. It
    cooperates with the FragmentDustMixDecorator class to handle each dust mix fragment separately
    (i.e. different grain materials and/or grain size bins). A subclass can specify a concrete
    destruction recipe by implementing a single function that takes properties of the radiation
    field and grain population/size bin under study and returns the corresponding non-destroyed
    density fraction.

    When a DustDestructionRecipe instance is included in the dynamic state configuration of a
    simulation, it requires there to be one or more medium components with a
    FragmentDustMixDecorator material mix that has the \em hasDynamicDensities flag enabled. These
    media components will be automatically detected and handled.

    <b>Convergence and noise</b>

    A cell counts as not converged if the density fraction of a fragment changes by more than the
    \em densityFractionTolerance. Near a destruction front, however, the Monte Carlo noise of the
    radiation field makes the fraction of some cells change in every iteration, by up to one for a
    sharp step, so that the number of such cells depends on the resolution of the front and on the
    number of photon packets rather than on the convergence of the solution. If the \em
    significance property is nonzero, the recipe therefore widens the tolerance of each cell by the
    change that the noise could cause: it calculates the density fraction for the radiation field
    multiplied by \f$1-k\sigma\f$ and by \f$1+k\sigma\f$, where \f$k\f$ is the significance and
    \f$\sigma\f$ the relative noise of the radiation field in the cell (see the
    MediumSystem::radiationFieldNoise() function), and adds the difference between these two
    fractions to the tolerance. This assumes that the fraction varies monotonically with the
    strength of the radiation field. A cell without any contributions to one of the noise batches
    has infinite noise, so that a change never counts against convergence. The medium system
    tracks the noise if this property is nonzero, at the cost of eight values per spatial cell. */
class DustDestructionRecipe : public DynamicStateRecipe
{
    ITEM_ABSTRACT(DustDestructionRecipe, DynamicStateRecipe, "a dust destruction recipe")

        PROPERTY_DOUBLE(densityFractionTolerance, "the convergence tolerance on the dynamic density fraction")
        ATTRIBUTE_MIN_VALUE(densityFractionTolerance, "[1e-6")
        ATTRIBUTE_MAX_VALUE(densityFractionTolerance, "0.5]")
        ATTRIBUTE_DEFAULT_VALUE(densityFractionTolerance, "0.05")

        PROPERTY_DOUBLE(significance, "the number of standard errors of the radiation field by which the "
                                      "convergence tolerance is widened (zero means no noise allowance)")
        ATTRIBUTE_MIN_VALUE(significance, "[0")
        ATTRIBUTE_MAX_VALUE(significance, "100]")
        ATTRIBUTE_DEFAULT_VALUE(significance, "0")
        ATTRIBUTE_DISPLAYED_IF(significance, "Level3")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function is called before an update cycle begins. It caches the radiation field
        wavelength grid and information on the dust population fragments in the medium system
        offered by one or more FragmentDustMixDecorator instances. */
    void beginUpdate(int numCells) override;

    /** This function is called repeatedly as part of the update cycle. If the density of the
        medium component in the cell under study is nonzero, it calculates the non-destroyed
        density fraction for each of the dust mix fragments in the medium component under study by
        calling the densityFraction() function defined in a subclass. If a newly calculated
        fraction differs by more than some small tolerance from the previously stored value, the
        corresponding medium state is updated. The function returns \em NotUpdated if none of the
        fractions have been updated, \em UpdatedNotConverged if one or more fractions were updated
        by more than the configured convergence tolerance (widened for the noise, if requested; see
        the class header), \em UpdatedConverged if one or more fractions were updated but none of
        them differed by more than that tolerance. */
    UpdateStatus update(MaterialState* state, const Array& Jv) override;

    /** This function returns true if the \em significance property is nonzero. */
    bool usesNoise() const override;

    /** This function returns the non-destroyed density fraction for a grain population with the
        specified type (graphite or silicate), average grain radius, radiation field, and
        equilibrium temperature. It must be implemented by a subclass. */
    virtual double densityFraction(bool graphite, double a, const Array& Jv, double T) const = 0;

    //======================== Private helpers =======================

private:
    /** This function returns the change of the density fraction of the specified fragment of the
        specified fragmented dust mix that the specified relative noise of the radiation field
        could cause, as described in the class header. */
    double noiseAllowance(const FragmentDustMixDecorator* fd, int f, const Array& Jv, double sigma) const;

    //======================== Data Members =======================

private:
    const MediumSystem* _ms{nullptr};              // the medium system, for the noise of the radiation field
    Array _dlambdav;                               // wavelength bin widths for the radiation field wavelength grid
    vector<const FragmentDustMixDecorator*> _fdv;  // fragmented mixes in the simulation (or nullptr for other mix)
};

////////////////////////////////////////////////////////////////////

#endif
