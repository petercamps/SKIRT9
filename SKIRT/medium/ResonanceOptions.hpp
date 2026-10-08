/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RESONANCEOPTIONS_HPP
#define RESONANCEOPTIONS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The ResonanceOptions class simply offers a number of configuration options related to the
    treatment of resonant line scattering, if this is enabled in the simulation. The options apply
    to all material mixes that support resonant line scattering, such as the
    LyaNeutralHydrogenGasMix class for the hydrogen Lyman-alpha line and the XRayIonicGasMix class
    for the resonance lines of hydrogen- and helium-like ions. */
class ResonanceOptions : public SimulationItem
{
    /** The enumeration type indicating the supported resonant line scattering acceleration schemes.

        In a high-optical-depth medium, the number of scatterings for photon packets near a
        resonance line such as Lyman-alpha is so high, and the corresponding free path lengths so
        short, that the Monte Carlo photon cycle becomes prohibitively slow. The core-skipping
        acceleration mechanism, first introduced for Lyman-alpha by Ahn et al. 2002 (ApJ,
        567:922-930) and used by many authors since, forces the wavelength of some photon packets
        from the core of the line into the wings of the line. This dramatically reduces the
        scattering cross section, allowing the photon packet to escape. More specifically, all
        photon packets with a dimensionless frequency \f$x\f$ smaller than a given critical value
        \f$x_\mathrm{crit}\f$ are treated this way.

        SKIRT implements three variations of the core-skipping acceleration scheme:

        - \em None: no acceleration is performed, corresponding to \f$x_\mathrm{crit}=0\f$. This
        can be useful for models with low optical depths, or to produce reference results for
        models with medium optical depths. For models with high optical depths, the run times will
        be prohibitively long.

        - \em Constant: acceleration with a constant critical value given by
        \f$x_\mathrm{crit}=3s\f$, where \f$s\f$ is the acceleration strength configured by the
        user. This option can be useful when the optical depth is fairly constant throughout the
        model, as is the case with many benchmark models.

        - \em Variable: acceleration with a variable critical value that depends on the local gas
        temperature and density. Specifically, the critical value is determined as \f[
        x_\mathrm{crit} = s\, \left( \frac{n_\mathrm{H}}{T} \right)^{1/6}, \f] where \f$s\f$ is the
        acceleration strength configured by the user, \f$n_\mathrm{H}\f$ is the number density (in
        \f$\mathrm{m}^{-3}\f$) of the medium component causing the scattering event (for
        Lyman-alpha, the neutral hydrogen number density), and \f$T\f$ is the gas temperature (in K)
        in the spatial cell hosting the scattering event. The rationale behind this formula is
        discussed below. This variable mechanism is applicable for most models, and is preferred
        for models with a broad dynamic range in optical depths.

        For both the constant and variable schemes, the user can configure the acceleration
        strength \f$s\f$, with a default value of unity. Larger values will decrease run time and
        accuracy; smaller values will increase run time and accuracy.

        <b>Rationale for the variable scheme</b>

        The approximate analytical solutions for the Lyman-alpha spectrum emerging from a static
        slab or sphere (Neufeld 1990, Dijkstra et al. 2006) depend on the product \f$a\tau_0\f$,
        where \f$a\f$ is the Voigt parameter and \f$\tau_0\f$ is the optical depth at the
        Lyman-alpha line center. Inspired by this result, many authors proposed acceleration
        schemes where \f$x_\mathrm{crit}\f$ is determined as a function of this product. For
        example, Smith et al. 2015 (MNRAS, 449, 4336-4362) used a critical value proportional to
        \f$(a\tau_0)^{1/3}\f$.

        However, calculating the optical depth requires selecting a path length. This has forced
        these schemes to depend on either a local scale (such as the size of the current spatial
        cell) or a global scale (such as the domain size). Both options seem undesirable, as they
        lead to a dependency on non-physical parameters (the resolution of the discretization or
        the portion of the physical world included in the model). Interestingly, Smith et al. 2015
        (MNRAS, 449, 4336-4362) noted that one could use the Jeans length as a physically motivated
        length scale. Expressing the Jeans length as well as the other quantities in
        \f$(a\tau_0)^{1/3}\f$ as a function of the local gas properties leads to \f$x_\mathrm{crit}
        \propto (n_\mathrm{H}/T)^{1/6}\f$. With the gas properties expressed in SI units,
        experiments with benchmark models show that a proportionality factor of order unity is
        appropriate. */
    ENUM_DEF(AccelerationScheme, None, Constant, Variable)
        ENUM_VAL(AccelerationScheme, None, "no acceleration")
        ENUM_VAL(AccelerationScheme, Constant, "acceleration scheme with a constant critical value")
        ENUM_VAL(AccelerationScheme, Variable, "acceleration scheme depending on local gas temperature and density")
    ENUM_END()

    ITEM_CONCRETE(ResonanceOptions, SimulationItem, "a set of options related to resonant line scattering")

        PROPERTY_ENUM(accelerationScheme, AccelerationScheme, "the resonant line scattering acceleration scheme")
        ATTRIBUTE_DEFAULT_VALUE(accelerationScheme, "Variable")
        ATTRIBUTE_DISPLAYED_IF(accelerationScheme, "Level2")

        PROPERTY_DOUBLE(accelerationStrength, "the acceleration strength; higher is faster but less accurate")
        ATTRIBUTE_MIN_VALUE(accelerationStrength, "]0")
        ATTRIBUTE_MAX_VALUE(accelerationStrength, "10]")
        ATTRIBUTE_DEFAULT_VALUE(accelerationStrength, "1")
        ATTRIBUTE_RELEVANT_IF(accelerationStrength, "accelerationSchemeConstant|accelerationSchemeVariable")
        ATTRIBUTE_DISPLAYED_IF(accelerationStrength, "Level2")

        PROPERTY_BOOL(includeHubbleFlow, "include the Doppler shift caused by the expansion of the universe")
        ATTRIBUTE_DEFAULT_VALUE(includeHubbleFlow, "false")
        ATTRIBUTE_DISPLAYED_IF(includeHubbleFlow, "Level2")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
