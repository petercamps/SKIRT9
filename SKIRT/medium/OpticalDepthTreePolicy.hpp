/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OPTICALDEPTHTREEPOLICY_HPP
#define OPTICALDEPTHTREEPOLICY_HPP

#include "MaterialWavelengthRangeInterface.hpp"
#include "TreePolicy.hpp"

//////////////////////////////////////////////////////////////////////

/** OpticalDepthTreePolicy is a tree subdivision policy that limits the optical depth across each
    cell for a given material type (dust, electrons, or gas), or for all media combined. A node is
    subdivided if the estimated extinction optical depth \f$\tau_\lambda\f$ at the configured
    wavelength \f$\lambda\f$ across the diagonal \f$\Delta s\f$ of the node exceeds the configured
    maximum optical depth \f$\tau_\text{max}\f$.

    Optical depth is additive. The policy therefore estimates the optical depth for each medium
    component of the selected material type (or for each medium component in the simulation for the
    \c All option), and adds these contributions: \f[ \tau_\lambda = \sum_h \kappa_{\lambda,h}\,
    \rho_h\, \Delta s, \f] where \f$\rho_h\f$ is the average density of medium component \f$h\f$
    inside the node, and \f$\kappa_{\lambda,h}\f$ is its extinction coefficient. For a dust
    component, \f$\rho_h\f$ is the mass density and \f$\kappa_{\lambda,h}\f$ the extinction cross
    section per unit mass, i.e. the cross section per hydrogen atom divided by the dust mass per
    hydrogen atom; for an electron or gas component, \f$\rho_h\f$ is the number density and
    \f$\kappa_{\lambda,h}\f$ the extinction cross section per entity. Using mass for dust allows
    sharing the density with the DensityTreePolicy. For the sake of performance, the extinction
    coefficient is assumed to be constant across the spatial domain. It is obtained from the default
    material mix of the medium component, i.e. the mix representative of the material properties of
    the component (see Medium::mix()), which is the configured mix for most media. The extinction
    cross section is evaluated for the default medium state, because the actual medium state is not
    yet known when the grid is constructed.

    The average density of a medium component inside a node is calculated exactly if the
    component offers the MassInBoxInterface, and is otherwise estimated from density samples (see
    the TreeNodeEvaluation class).

    If the simulation has no medium components of the selected material type, or no medium
    components at all for the \c All option, setup reports a fatal error.

    This class implements the MaterialWavelengthRangeInterface to indicate that
    wavelength-dependent material properties will be required for the configured wavelength. */
class OpticalDepthTreePolicy : public TreePolicy, public MaterialWavelengthRangeInterface
{
    /** The enumeration type indicating the material type to which the criterion applies. */
    ENUM_DEF(MaterialType, Dust, Electrons, Gas, All)
        ENUM_VAL(MaterialType, Dust, "dust")
        ENUM_VAL(MaterialType, Electrons, "electrons")
        ENUM_VAL(MaterialType, Gas, "gas")
        ENUM_VAL(MaterialType, All, "all media combined")
    ENUM_END()

    ITEM_CONCRETE(OpticalDepthTreePolicy, TreePolicy,
                  "a tree subdivision policy limiting the optical depth across a cell")
        ATTRIBUTE_TYPE_DISPLAYED_IF(OpticalDepthTreePolicy, "Level2")

        PROPERTY_ENUM(materialType, MaterialType, "the material type to which the criterion applies")
        ATTRIBUTE_DEFAULT_VALUE(materialType, "DustMix:Dust;All")

        PROPERTY_DOUBLE(maxOpticalDepth, "the maximum diagonal optical depth for each cell")
        ATTRIBUTE_MIN_VALUE(maxOpticalDepth, "]0")
        ATTRIBUTE_MAX_VALUE(maxOpticalDepth, "100]")

        PROPERTY_DOUBLE(wavelength, "the wavelength at which to evaluate the optical depth")
        ATTRIBUTE_QUANTITY(wavelength, "wavelength")
        ATTRIBUTE_MIN_VALUE(wavelength, "1 pm")
        ATTRIBUTE_MAX_VALUE(wavelength, "1 m")
        ATTRIBUTE_DEFAULT_VALUE(wavelength, "0.55 micron")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function obtains and caches information used by the needsSubdivide() function. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns true if the estimated optical depth across the diagonal of the
        specified node exceeds the configured maximum optical depth, and false otherwise. */
    bool needsSubdivide(TreeNodeEvaluation& node) const override;

    /** This function returns a wavelength range corresponding to the configured wavelength,
        indicating that wavelength-dependent material properties will be required for this
        wavelength. */
    Range wavelengthRange() const override;

    //======================== Data Members ========================

private:
    // data members initialized by setupSelfBefore()
    vector<int> _hv;         // the indices of the medium components included in the criterion
    vector<double> _kappav;  // the extinction coefficient for each of these medium components
};

//////////////////////////////////////////////////////////////////////

#endif
