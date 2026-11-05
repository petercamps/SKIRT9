/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TEMPERATUREWAVELENGTHCELLLIBRARY_HPP
#define TEMPERATUREWAVELENGTHCELLLIBRARY_HPP

#include "SpatialCellLibrary.hpp"

//////////////////////////////////////////////////////////////////////

/** The TemperatureWavelengthCellLibrary class provides a library scheme for grouping spatial cells
    based on an indicative dust temperature and an indicative dust wavelength of the stored
    radiation field in each cell. The library contains a two-dimensional table of entries
    corresponding to different values for each of these properties. The user can configure the
    number of bins in each of the two dimensions.

    The cells mapped to a library entry share the mean intensity of the radiation field averaged
    over these cells, while the emission of each cell is still calculated for the dust in that
    cell. The properties used for the mapping must therefore characterize the radiation field,
    and they must not depend on the dust composition. With multiple dust medium components that
    have different spatial distributions, the composition varies from cell to cell. If the
    properties depended on it, cells with different radiation fields would share an entry, and
    the averaged radiation field would misrepresent the field in each of them, an error that
    more library entries do not remove. Within a single component, the composition cannot vary,
    because a cell library cannot be combined with a material mix family or with a mix that has
    extra specific state variables, such as a fragmented dust mix.

    To this end, both properties are calculated for the dust in the simulation as a whole rather
    than for the dust in the cell, weighting each dust medium component \f$h\f$ by its total
    mass \f$M_h\f$ in the simulation. The indicative dust temperature for a cell \f$m\f$ is the
    mass-weighted average \f[ {\bar{T}}_m = \frac{\sum_h M_h\, T_h(J_m)}{\sum_h M_h} \f] where
    \f$T_h(J_m)\f$ is the equilibrium temperature of the representative grain of component
    \f$h\f$ embedded in the radiation field \f$J_m\f$ of the cell. The indicative dust
    wavelength for the cell is defined as \f[ {\bar{\lambda}}_m = \frac{\int_0^\infty
    \bar\kappa_\lambda^{\text{abs}}\, J_{m,\lambda}\, \lambda\, {\text{d}}\lambda }{
    \int_0^\infty \bar\kappa_\lambda^{\text{abs}}\, J_{m,\lambda}\, {\text{d}}\lambda} \f] where
    \f[ \bar\kappa_\lambda^{\text{abs}} = \frac{\sum_h M_h\,
    \kappa_{h,\lambda}^{\text{abs}}}{\sum_h M_h} \f] is the mass-weighted average of the
    absorption mass coefficients of the dust components. Both properties thus depend on the
    radiation field of the cell only. If the total dust mass in the simulation is zero, all
    components receive the same weight.

    If the dust has the same composition in every cell, for example with a single dust
    component, these properties equal the indicative dust temperature of the dust in the cell
    (see the MediumSystem::indicativeDustTemperature() function) and its mean wavelength
    weighted by the absorption opacity of that dust. */
class TemperatureWavelengthCellLibrary : public SpatialCellLibrary
{
    ITEM_CONCRETE(TemperatureWavelengthCellLibrary, SpatialCellLibrary,
                  "a library scheme for grouping spatial cells based on indicative temperature and wavelength")
        ATTRIBUTE_TYPE_INSERT(TemperatureWavelengthCellLibrary, "NonIdentitySpatialCellLibrary")

        PROPERTY_INT(numTemperatures, "the number of temperature bins")
        ATTRIBUTE_MIN_VALUE(numTemperatures, "5")
        ATTRIBUTE_MAX_VALUE(numTemperatures, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numTemperatures, "40")

        PROPERTY_INT(numWavelengths, "the number of wavelength bins")
        ATTRIBUTE_MIN_VALUE(numWavelengths, "5")
        ATTRIBUTE_MAX_VALUE(numWavelengths, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numWavelengths, "25")

    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function returns the number of entries in the library. In this class the function
        returns the user-configured number of bins, i.e. the product of the number of temperature
        bins and the number of wavelength bins. */
    int numEntries() const override;

    /** This function returns a vector \em nv with length \f$N_{\text{cells}}\f$ that maps each
        cell index \f$m\f$ to the corresponding library entry index \f$n_m\f$. In this class the
        function loops over all spatial cells and calculates the indicative dust temperature and
        the indicative dust wavelength of the stored radiation field for each of them, as
        described in the class header. Based on
        these values, a two-dimensional grid is established such that it fits all the measured
        values. The temperature grid points \f${\bar{T}}_{(i)}\f$ are distributed linearly, i.e.
        \f[ {\bar{T}}_{(i)} = {\bar{T}}_{\text{min}} + \frac{i}{N_{\bar{T}}}\,
        ({\bar{T}}_{\text{max}} - {\bar{T}}_{\text{min}}) \qquad i=0,\ldots,N_{\bar{T}} \f] where
        \f${\bar{T}}_{\text{min}}\f$ and \f${\bar{T}}_{\text{max}}\f$ represent the smallest and
        largest values of the indicative dust temperature found among all spatial cells. The
        wavelength grid points have logarithmic distribution, \f[ {\bar{\lambda}}_{(j)} =
        {\bar{\lambda}}_{\text{min}} \left( \frac{ {\bar{\lambda}}_{\text{max}} }{
        {\bar{\lambda}}_{\text{min}} } \right)^{j/N_{\bar{\lambda}}} \qquad
        j=0,\ldots,N_{\bar{\lambda}}. \f] The function then calculates for each cell \f$m\f$ its
        library entry \f$n \equiv (i,j)\f$. */
    vector<int> mapping(const Array& bv) const override;
};

////////////////////////////////////////////////////////////////////

#endif
