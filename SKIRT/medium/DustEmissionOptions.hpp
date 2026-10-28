/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef DUSTEMISSIONOPTIONS_HPP
#define DUSTEMISSIONOPTIONS_HPP

#include "DisjointWavelengthGrid.hpp"
#include "SimulationItem.hpp"
#include "SourceWavelengthRangeInterface.hpp"
#include "SpatialCellLibrary.hpp"
#include "WavelengthDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** The DustEmissionOptions class simply offers a number of configuration options related to
    thermal emission from dust. In a mode where dust emission is enabled, the simulation also
    needs a wavelength grid on which to calculate the dust emission spectrum.

    <b>Convergence of secondary emission iterations</b>

    If the simulation iterates over secondary emission, the dust absorbs part of its own emission,
    and is heated by it. Several options control when the iteration loop is considered to have
    converged with regard to this dust self-heating. In each iteration, the simulation determines
    the absorbed primary luminosity \f$L_\mathrm{p}\f$, i.e. the luminosity of the primary
    radiation absorbed by the dust; the absorbed secondary luminosity \f$L_\mathrm{s}(n)\f$, i.e.
    the luminosity of the dust emission launched in iteration \f$n\f$ that is absorbed by the dust;
    and the emitted dust luminosity \f$L_\mathrm{e}(n)\f$. It then calculates the luminosity
    deficit \f$\delta(n) = (L_\mathrm{s}(n) - L_\mathrm{e}(n) + L_\mathrm{p}) / L_\mathrm{p}\f$,
    i.e. the relative amount by which the escaping dust luminosity falls short of the absorbed
    primary luminosity, which it must equal once the solution has converged, and the secondary
    change \f$c(n) = |L_\mathrm{s}(n) - L_\mathrm{s}(n-k)| / (k\,L_\mathrm{s}(n))\f$ with
    \f$k=\min(w,n-2)\f$, i.e. the mean relative change of the absorbed secondary luminosity per
    iteration over the last \f$w\f$ iterations, where \f$w\f$ is the value of the \em
    numConvergenceIterations property, excluding the first iteration of the loop, whose dust
    emission lacks any self-heating.

    The iteration over dust emission has converged if the absorbed secondary luminosity is smaller
    than the fraction \em maxFractionOfPrimary of the absorbed primary luminosity, or if at least
    \f$w+1\f$ iterations and at least 3 iterations have been performed, the absolute value of the
    mean luminosity deficit over the last \f$w\f$ iterations is at most \em maxLuminosityDeficit,
    and the secondary change is at most \em maxSecondaryChange. The luminosity deficit is
    meaningful only if thermal dust emission is the only form of secondary emission. If the
    simulation also has gas emission, the condition on the deficit is omitted.

    <b>Rebalance of the secondary radiation field</b>

    If the \em rebalanceSecondaryEmission property is true (the default), the simulation iterates
    over secondary emission, and thermal dust emission is the only form of secondary emission, the
    secondary radiation field is rebalanced after each iteration from the second one on. The
    spatial cells are grouped into regions based on the ratio of the absorbed secondary to the
    absorbed primary luminosity and on the absorbed luminosity per unit dust mass, with the number
    of bins on each axis given by the \em numAbsorptionRatioBins and \em numHeatingRateBins
    properties. The simulation records the luminosity transfer between the regions during the
    iteration, solves for the dust emission of each region that is consistent with this transfer,
    and rescales the secondary radiation field in each region accordingly. This greatly reduces the
    number of iterations needed for models with strong dust self-heating. For more information, see
    the SelfHeatingRebalancer class. */
class DustEmissionOptions : public SimulationItem, public SourceWavelengthRangeInterface
{
    /** The enumeration type indicating the method used for dust emission calculations. */
    ENUM_DEF(EmissionType, Equilibrium, Stochastic)
        ENUM_VAL(EmissionType, Equilibrium, "assume local thermal equilibrium (LTE)")
        ENUM_VAL(EmissionType, Stochastic, "handle stochastically heated grains (non-LTE)")
    ENUM_END()

    ITEM_CONCRETE(DustEmissionOptions, SimulationItem, "a set of options related to thermal emission from dust")

        PROPERTY_ENUM(dustEmissionType, EmissionType, "the method used for dust emission calculations")
        ATTRIBUTE_DEFAULT_VALUE(dustEmissionType, "Equilibrium")
        ATTRIBUTE_INSERT(dustEmissionType, "dustEmissionTypeStochastic:StochasticDustEmission")
        ATTRIBUTE_DISPLAYED_IF(dustEmissionType, "Level2")

        PROPERTY_BOOL(includeHeatingByCMB, "add the cosmic microwave background (CMB) as a dust heating source term")
        ATTRIBUTE_DEFAULT_VALUE(includeHeatingByCMB, "false")
        ATTRIBUTE_DISPLAYED_IF(includeHeatingByCMB, "NonZeroRedshift")

        PROPERTY_ITEM(cellLibrary, SpatialCellLibrary, "the spatial cell grouping scheme for calculating dust emission")
        ATTRIBUTE_DEFAULT_VALUE(cellLibrary, "AllCellsLibrary")
        ATTRIBUTE_REQUIRED_IF(cellLibrary, "false")
        ATTRIBUTE_DISPLAYED_IF(cellLibrary, "Level2")

        PROPERTY_ITEM(dustEmissionWLG, DisjointWavelengthGrid,
                      "the wavelength grid for calculating the dust emission spectrum")
        ATTRIBUTE_DEFAULT_VALUE(dustEmissionWLG, "LogWavelengthGrid")

        PROPERTY_DOUBLE(maxFractionOfPrimary, "convergence is reached when the total absorbed dust luminosity "
                                              "is less than this fraction of the total absorbed primary luminosity")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPrimary, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPrimary, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPrimary, "0.01")
        ATTRIBUTE_RELEVANT_IF(maxFractionOfPrimary, "IterateSecondary")

        PROPERTY_DOUBLE(maxLuminosityDeficit, "convergence requires the mean deficit of the escaping dust "
                                              "luminosity relative to the absorbed primary luminosity to be smaller")
        ATTRIBUTE_MIN_VALUE(maxLuminosityDeficit, "]0")
        ATTRIBUTE_MAX_VALUE(maxLuminosityDeficit, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxLuminosityDeficit, "0.01")
        ATTRIBUTE_RELEVANT_IF(maxLuminosityDeficit, "IterateSecondary")

        PROPERTY_DOUBLE(maxSecondaryChange, "convergence requires the mean relative change per iteration "
                                            "of the absorbed secondary luminosity to be smaller")
        ATTRIBUTE_MIN_VALUE(maxSecondaryChange, "]0")
        ATTRIBUTE_MAX_VALUE(maxSecondaryChange, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxSecondaryChange, "0.01")
        ATTRIBUTE_RELEVANT_IF(maxSecondaryChange, "IterateSecondary")

        PROPERTY_INT(numConvergenceIterations, "the number of iterations over which convergence is evaluated")
        ATTRIBUTE_MIN_VALUE(numConvergenceIterations, "1")
        ATTRIBUTE_MAX_VALUE(numConvergenceIterations, "10")
        ATTRIBUTE_DEFAULT_VALUE(numConvergenceIterations, "3")
        ATTRIBUTE_RELEVANT_IF(numConvergenceIterations, "IterateSecondary")

        PROPERTY_BOOL(rebalanceSecondaryEmission,
                      "rebalance the secondary radiation field over regions after each iteration")
        ATTRIBUTE_DEFAULT_VALUE(rebalanceSecondaryEmission, "true")
        ATTRIBUTE_RELEVANT_IF(rebalanceSecondaryEmission, "IterateSecondary")

        PROPERTY_INT(numAbsorptionRatioBins,
                     "the number of rebalance bins in the ratio of secondary to primary absorption")
        ATTRIBUTE_MIN_VALUE(numAbsorptionRatioBins, "1")
        ATTRIBUTE_MAX_VALUE(numAbsorptionRatioBins, "100")
        ATTRIBUTE_DEFAULT_VALUE(numAbsorptionRatioBins, "10")
        ATTRIBUTE_RELEVANT_IF(numAbsorptionRatioBins, "IterateSecondary&rebalanceSecondaryEmission")
        ATTRIBUTE_DISPLAYED_IF(numAbsorptionRatioBins, "Level3")

        PROPERTY_INT(numHeatingRateBins, "the number of rebalance bins in the absorbed luminosity per unit dust mass")
        ATTRIBUTE_MIN_VALUE(numHeatingRateBins, "1")
        ATTRIBUTE_MAX_VALUE(numHeatingRateBins, "100")
        ATTRIBUTE_DEFAULT_VALUE(numHeatingRateBins, "10")
        ATTRIBUTE_RELEVANT_IF(numHeatingRateBins, "IterateSecondary&rebalanceSecondaryEmission")
        ATTRIBUTE_DISPLAYED_IF(numHeatingRateBins, "Level3")

        PROPERTY_DOUBLE(sourceWeight, "the weight of dust emission for the number of photon packets launched")
        ATTRIBUTE_MIN_VALUE(sourceWeight, "]0")
        ATTRIBUTE_MAX_VALUE(sourceWeight, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(sourceWeight, "1")
        ATTRIBUTE_DISPLAYED_IF(sourceWeight, "Level3")

        PROPERTY_DOUBLE(wavelengthBias,
                        "the fraction of secondary photon packet wavelengths sampled from a bias distribution")
        ATTRIBUTE_MIN_VALUE(wavelengthBias, "[0")
        ATTRIBUTE_MAX_VALUE(wavelengthBias, "1]")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthBias, "0.5")
        ATTRIBUTE_DISPLAYED_IF(wavelengthBias, "Level3")

        PROPERTY_ITEM(wavelengthBiasDistribution, WavelengthDistribution,
                      "the bias distribution for sampling secondary photon packet wavelengths")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthBiasDistribution, "DefaultWavelengthDistribution")
        ATTRIBUTE_RELEVANT_IF(wavelengthBiasDistribution, "wavelengthBias")
        ATTRIBUTE_DISPLAYED_IF(wavelengthBiasDistribution, "Level3")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the wavelength range of the dust emission wavelength grid,
        implementing the SourceWavelengthRangeInterface interface so that the wavelength bias
        distribution configured here can locate its associated "source" wavelength range. */
    Range wavelengthRange() const override { return _dustEmissionWLG->wavelengthRange(); }
};

////////////////////////////////////////////////////////////////////

#endif
