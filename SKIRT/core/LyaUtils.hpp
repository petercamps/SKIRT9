/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef LYAUTILS_HPP
#define LYAUTILS_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** This namespace offers utility functions related to the treatment of Lyman-alpha (Lya) line
    transfer.

    When a Lya photon is absorbed by a neutral hydrogen atom in the ground state, the atom is
    excited to the 2p energy level and a new Lya photon is emitted almost immediately as a result
    of the subsequent downward transition. This happens fast enough that we can consider the
    combined process as a scattering event.

    The cross section for a single hydrogen atom of the Lya scattering of a photon can be derived
    using quantum mechanical considerations, resulting in a sharply peaked profile as a function of
    the photon wavelength in the atom's rest frame. Because each atom has its own velocity, a
    photon with a given wavelength in the gas rest frame will appear Doppler shifted to a slightly
    different wavelength for each atom in the gas. To compute the Lya absorption cross section for
    a collection of moving atoms, we must therefore convolve the single-atom cross section with the
    atom velocity distribution, which in turn depends on the gas temperature.

    Assuming a Maxwell-Boltzmann velocity distribution, we define the characteristic thermal
    velocity \f$v_\mathrm{th}\f$ as \f[ v_\mathrm{th} = \sqrt{\frac{2 k_\mathrm{B}
    T}{m_\mathrm{p}}} \f] where \f$k_\mathrm{B}\f$ is the Boltzmann constant, \f$m_\mathrm{p}\f$ is
    the proton mass, and \f$T\f$ is the temperature of the gas. We then introduce the dimensionless
    frequency variable \f$x\f$, defined as \f[ x = \frac{\nu - \nu_\alpha}{\nu_\alpha}
    \,\frac{c}{v_\mathrm{th}} = \frac{v_\mathrm{p}}{v_\mathrm{th}} \f]

    where \f$\nu=c/\lambda\f$ is the regular frequency variable, \f$\nu_\alpha=c/\lambda_\alpha\f$
    is the frequency at the Lya line center, \f$\lambda_\alpha\f$ is the wavelength at the Lya line
    center, and \f$c\f$ is the speed of light in vacuum. The last equality introduces the velocity
    shift \f$v_\mathrm{p}\f$ of the photon frequency relative to the Lya line center, defined by
    \f[ \frac{v_\mathrm{p}}{c} = \frac{\nu - \nu_\alpha}{\nu_\alpha} \approx -\frac{\lambda -
    \lambda_\alpha}{\lambda_\alpha} \f] where the approximate equality holds for \f$v_\mathrm{p}\ll
    c\f$.

    After neglecting some higher order terms, the convolution of the single-atom profile with the
    Maxwell-Boltzmann velocity distribution yields the following expression for the
    velocity-weighted Lya scattering cross section \f$\sigma_\alpha(x,T)\f$ of a hydrogen atom in
    gas at temperature \f$T\f$ as a function of the dimensionless photon frequency \f$x\f$: \f[
    \sigma_\alpha(x,T) = \sigma_{\alpha,0}(T)\,H(a_\mathrm{v}(T),x) \f] where the cross section at
    the line center \f$\sigma_{\alpha,0}(T)\f$ is given by \f[ \sigma_{\alpha,0}(T) =
    \frac{3\lambda_\alpha^2 a_\mathrm{v}(T)}{2\sqrt{\pi}}, \f] the Voigt parameter
    \f$a_\mathrm{v}(T)\f$ is given by \f[ a_\mathrm{v}(T) = \frac{A_\alpha}{4\pi\nu_\alpha}
    \,\frac{c}{v_\mathrm{th}} \f] with \f$A_\alpha\f$ the Einstein A-coefficient of the Lya
    transition; and the Voigt function \f$H(a_\mathrm{v},x)\f$ is defined by \f[ H(a_\mathrm{v},x)
    = \frac{a_\mathrm{v}}{\pi} \int_{-\infty}^\infty \frac{\mathrm{e}^{-y^2} \,\mathrm{d}y}
    {(y-x)^2+a_\mathrm{v}^2} \f] which is normalized so that \f$H(a_\mathrm{v},0) = 1\f$. */
namespace LyaUtils
{
    /** This function returns the Lyman-alpha scattering cross section per hydrogen atom
        \f$\sigma_\alpha(x, T)\f$ at the given dimensionless photon frequency and gas temperature.
        */
    double sectionForDimlessFreq(double x, double T);

    /** This function returns the Lyman-alpha scattering cross section per hydrogen atom
        \f$\sigma_\alpha(\lambda, T)\f$ at the given photon wavelength and gas temperature. */
    double sectionForWavelength(double lambda, double T);
}

////////////////////////////////////////////////////////////////////

#endif