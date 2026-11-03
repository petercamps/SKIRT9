/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RADIATIONFIELDGRADIENTCRITERION_HPP
#define RADIATIONFIELDGRADIENTCRITERION_HPP

#include "GradientCriterion.hpp"
class DisjointWavelengthGrid;
class MediumSystem;

////////////////////////////////////////////////////////////////////

/** A RadiationFieldGradientCriterion flags the spatial cells in which the mean intensity of the
    radiation field, integrated over a given wavelength range, differs strongly from that in a
    neighboring cell, as described for the GradientCriterion class. It is used for dynamic grid
    refinement (see the DynamicRefinementOptions class), for example to resolve the attenuation of
    the radiation entering an optically thick region. Unlike the DustTemperatureGradientCriterion,
    it does not involve the properties of the material, and its field is cheap to calculate.

    The field of a cell is the integrated mean intensity \f$J = \sum_\ell (J_\lambda)_\ell \,
    (\Delta\lambda)_\ell\f$, summed over the bins \f$\ell\f$ of the radiation field wavelength grid
    whose characteristic wavelength lies in the range from \em minWavelength to \em maxWavelength
    (see the MediumSystem::meanIntensity() function). By default, the range includes all bins. The
    field is expressed in the internal units of SKIRT for this quantity, i.e. \f$\text{W}\,
    \text{m}^{-2}\,\text{sr}^{-1}\f$. Because the mean intensity varies by orders of magnitude, the
    criterion is normally used with the \em logarithmic flag, so that the threshold is a relative
    change. A temperature in equilibrium with the radiation field varies roughly as the fourth to
    sixth root of the absorbed intensity, so that a relative change of the mean intensity
    corresponds to a four to six times smaller relative change of such a temperature, to the extent
    that the spectral shape of the radiation field does not change between the cells. Setup reports
    a fatal error if the simulation does not store the radiation field, or if the wavelength range
    contains no bin of the radiation field wavelength grid.

    The noise of the field, used if the \em significance property is nonzero (see the
    GradientCriterion class), is estimated as the field value multiplied by the relative noise of
    the radiation field of the cell summed over all wavelength bins (see the
    MediumSystem::radiationFieldNoise() function). For a narrow integration range, this
    underestimates the noise, because fewer photon packets contribute to it.

    The criterion takes part in all iteration loops, because each loop recalculates at least part of
    the radiation field. In the secondary emission loop, however, the primary radiation field is not
    recalculated, as described for the DustTemperatureGradientCriterion class. The number density
    for the column density limit is the sum of the number densities of all medium components. */
class RadiationFieldGradientCriterion : public GradientCriterion
{
    ITEM_CONCRETE(RadiationFieldGradientCriterion, GradientCriterion,
                  "a criterion on the difference of the integrated mean intensity between neighboring cells")
        ATTRIBUTE_TYPE_ALLOWED_IF(RadiationFieldGradientCriterion, "RadiationField")

        PROPERTY_DOUBLE(minWavelength, "the shortest wavelength of the integration range")
        ATTRIBUTE_QUANTITY(minWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(minWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(minWavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(minWavelength, "1 pm")

        PROPERTY_DOUBLE(maxWavelength, "the longest wavelength of the integration range")
        ATTRIBUTE_QUANTITY(maxWavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(maxWavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(maxWavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(maxWavelength, "1 m")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the simulation stores the radiation field, and determines the
        bins of the radiation field wavelength grid in the integration range. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns true for any iteration loop. */
    bool changesIn(IterationHistory::Loop loop) const override;

    /** This function returns the integrated mean intensity in the spatial cell with index \f$m\f$.
        */
    double value(int m) const override;

    /** This function returns an estimate of the noise of the integrated mean intensity in the
        spatial cell with index \f$m\f$, as described in the class header. */
    double noise(int m) const override;

    /** This function returns the sum of the number densities of all medium components in the
        spatial cell with index \f$m\f$. */
    double numberDensity(int m) const override;

    /** This function returns a description of the field. */
    string fieldDescription() const override;

    //======================== Data Members ========================

private:
    const MediumSystem* _ms{nullptr};
    const DisjointWavelengthGrid* _wavelengthGrid{nullptr};
    vector<int> _bins;  // indices of the radiation field wavelength bins in the integration range
};

////////////////////////////////////////////////////////////////////

#endif
