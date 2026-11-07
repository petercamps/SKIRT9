/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ESCAPEINSTRUMENT_HPP
#define ESCAPEINSTRUMENT_HPP

#include "EscapeBinning.hpp"
#include "Instrument.hpp"

////////////////////////////////////////////////////////////////////

/** An EscapeInstrument records the radiation that escapes from the model, in bins of escape
    direction and wavelength. Rather than receiving peel-off photon packets, it records at each
    step of each photon packet the part of the packet's luminosity that leaves the model along the
    packet's own direction without further interaction. With forced scattering, this is the
    fraction \f$\exp(-\tau_\mathrm{path})\f$ that the simulation removes from the packet's weight at
    each step, where \f$\tau_\mathrm{path}\f$ is the optical depth of the complete path to the model
    boundary. Without forced scattering, it is the luminosity of a photon packet whose interaction
    point lies beyond the model boundary, attenuated by the absorption optical depth of the
    remaining path if absorption is treated explicitly. In a simulation without media, every photon
    packet escapes unattenuated at launch. As for peel-off instruments, recording happens only
    during the peel-off segments of the simulation, for primary and secondary emission alike.

    The recorded luminosity is an unbiased estimate of the radiation escaping in each bin, which
    requires no extra path calculation. For wide bins, the noise is low because every step of every
    photon packet contributes to some bin. For individual lines of sight, peel-off instruments are
    more efficient. The two kinds of instruments can be combined in a single simulation.

    The binning of the escape directions is configured through an EscapeBinning: all directions
    together, bins of inclination and azimuth, or HEALPix pixels. For a bin that covers a fraction
    \f$f\f$ of the sky, the instrument reports the flux density \f$F_\lambda =
    L_\mathrm{bin} / (4\pi d^2 f \Delta\lambda)\f$, the average flux density received by
    observers at distance \f$d\f$ in the directions of the bin. As for distant instruments, a zero
    distance places the instrument in the observer frame corresponding to the model redshift (see
    DistantInstrument).

    The instrument records the same flux components and statistics as an SED instrument, except
    for polarization, which is not supported. With the \em SEDFiles output style, the instrument
    writes the spectrum of each bin to a separate SED column file; with the \em DataCubes output
    style, it writes a FITS data cube for each flux component, with the bins along the spatial axes
    and the wavelengths along the third axis. With either style, it also writes a column text file
    <tt>prefix_instr_bins.dat</tt> listing the index, the sky fraction, and the defining values of
    each bin. */
class EscapeInstrument : public Instrument
{
    /** The enumeration type indicating the output style of an escape instrument. */
    ENUM_DEF(OutputStyle, SEDFiles, DataCubes)
        ENUM_VAL(OutputStyle, SEDFiles, "an SED file for each escape bin")
        ENUM_VAL(OutputStyle, DataCubes, "a data cube for each flux component, with the bins along the spatial axes")
    ENUM_END()

    ITEM_CONCRETE(EscapeInstrument, Instrument, "an escape instrument (records the radiation escaping from the model)")
        ATTRIBUTE_TYPE_DISPLAYED_IF(EscapeInstrument, "Level2")

        PROPERTY_DOUBLE(distance, "the distance to the system")
        ATTRIBUTE_QUANTITY(distance, "distance")
        ATTRIBUTE_MIN_VALUE(distance, "[0")

        PROPERTY_ITEM(binning, EscapeBinning, "the binning of the escape directions")
        ATTRIBUTE_DEFAULT_VALUE(binning, "AllDirectionsEscapeBinning")

        PROPERTY_ENUM(outputStyle, OutputStyle, "the output style")
        ATTRIBUTE_DEFAULT_VALUE(outputStyle, "SEDFiles")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function configures the flux recorder with the appropriate frame and distances. */
    void setupSelfBefore() override;

    /** This function configures the flux recorder for the escape bins provided by the binning,
        which has been set up at this point, before the base class finalizes the configuration. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function records the luminosity of the specified photon packet that escapes from the
        model along the packet's direction, i.e. the packet's luminosity multiplied by the
        specified escape fraction, in the bin of the packet's direction. If the \em
        includeTransparent flag is false, the transparent contribution of a photon packet that has
        not been scattered is not recorded (see FluxRecorder::detect()). The function does nothing
        if the instrument does not record the packet's wavelength or direction. It is called from
        many threads at once and does not change any state of its own. */
    void detectEscape(const PhotonPacket* pp, double escapeFraction, bool includeTransparent);

    /** This function calibrates the instrument and outputs the recorded spectra, as described for
        the base class, and writes the bin table. */
    void write() override;
};

////////////////////////////////////////////////////////////////////

#endif
