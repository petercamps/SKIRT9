/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INSTRUMENT_HPP
#define INSTRUMENT_HPP

#include "SimulationItem.hpp"
#include "WavelengthGrid.hpp"
class FluxRecorder;
class PhotonPacket;

////////////////////////////////////////////////////////////////////

/** Instrument is an abstract class representing instruments to collect the photon packets during a
    Monte Carlo simulation. There are two kinds of instruments. Peel-off instruments (subclasses of
    PeelOffInstrument) record a peel-off photon packet sent towards them at each emission and
    scattering event. Escape instruments (EscapeInstrument) record the luminosity that escapes
    from the model along the photon packet's own direction at each step. The simulation serves the
    two kinds of instruments separately.

    This top-level abstract class offers the facilities that are common to all instruments: the
    instrument name, the wavelength grid, the user properties that configure which flux
    contributions need to be recorded, and the FluxRecorder instance that actually records the
    detected photon packets. A wavelength grid is established either by specifying a grid for this
    instrument specifically, or by defaulting to the default grid of the wavelength grid pool. */
class Instrument : public SimulationItem
{
    ITEM_ABSTRACT(Instrument, SimulationItem, "an instrument")

        PROPERTY_STRING(instrumentName, "the name for this instrument")

        PROPERTY_ITEM(wavelengthGrid, WavelengthGrid, "the wavelength grid for this instrument")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthGrid, "!DefaultInstrumentWavelengthGrid:LogWavelengthGrid;")
        ATTRIBUTE_RELEVANT_IF(wavelengthGrid, "Panchromatic")
        ATTRIBUTE_REQUIRED_IF(wavelengthGrid, "!DefaultInstrumentWavelengthGrid")
        ATTRIBUTE_DISPLAYED_IF(wavelengthGrid, "Level2|!DefaultInstrumentWavelengthGrid")

        ATTRIBUTE_SUB_PROPERTIES_HERE(Instrument)

        PROPERTY_BOOL(recordComponents, "record flux components separately")
        ATTRIBUTE_DEFAULT_VALUE(recordComponents, "false")

        PROPERTY_INT(numScatteringLevels, "the number of individually recorded scattering levels")
        ATTRIBUTE_MIN_VALUE(numScatteringLevels, "0")
        ATTRIBUTE_MAX_VALUE(numScatteringLevels, "99")
        ATTRIBUTE_DEFAULT_VALUE(numScatteringLevels, "0")
        ATTRIBUTE_RELEVANT_IF(numScatteringLevels, "recordComponents")
        ATTRIBUTE_DISPLAYED_IF(numScatteringLevels, "Level2")

        PROPERTY_BOOL(recordStatistics, "record information for calculating statistical properties")
        ATTRIBUTE_DEFAULT_VALUE(recordStatistics, "false")
        ATTRIBUTE_DISPLAYED_IF(recordStatistics, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function performs setup for the instrument. It establishes the wavelength grid for the
        instrument: if a grid is specified for the instrument, that grid is used (or the referenced
        grid in the case of a ReferenceWavelengthGrid). If not, the default wavelength grid of the
        WavelengthGridPool is used instead. If neither of these grids are specified, the function
        throws a fatal error.

        The function also creates and partially configures the FluxRecorder instance for this
        instrument, passing it the values of the user properties offered by this class and some
        extra information on the simulation. The setup functions of each subclass are expected to
        augment the configuration by calling the appropriate configuration functions of the flux
        recorder. */
    void setupSelfBefore() override;

    /** This function finalizes the configuration of FluxRecorder instance for this instrument. */
    void setupSelfAfter() override;

    /** The destructor releases the FluxRecorder instance for this instrument. */
    ~Instrument();

    /** This function configures the flux recorder of a distant instrument with the specified
        distance, intended for use in the setup of subclasses. If the distance is nonzero, the
        instrument is placed in the rest frame of the model at that (non-relativistic) distance. If
        the distance is zero, the instrument is placed in an observer frame corresponding to the
        redshift of the model, using the angular-diameter and luminosity distances derived from the
        cosmology. If the distance and the model redshift are both zero, the function throws a fatal
        error. */
    void configureRecorderDistance(double distance);

    //======================== Other Functions =======================

public:
    /** This function returns the instrument name as human-readable name for the simulation item,
        so that it can be used in log messages to identify the instrument and differentiate it from
        other instruments. */
    string itemName() const override;

    /** This function returns the wavelength grid for the instrument as determined during setup,
        i.e. either the grid specified for this instrument (or the grid it references) or the
        default grid of the wavelength grid pool. After setup has completed, the function never
        returns a nulltpr because setupSelfBefore() throws a fatal error if neither of these grids
        are specified. */
    const WavelengthGrid* instrumentWavelengthGrid() const { return _instrumentWavelengthGrid; }

    /** This function flushes any information buffered by the flux recorder during detection. It
        simply calls the corresponding function of the FluxRecorder instance associated with this
        instrument. */
    void flush();

    /** This function calibrates the instrument and outputs the recorded contents to a set of
        files. The implementation in this class calls the corresponding function of the
        FluxRecorder instance associated with this instrument. Subclasses may override this function
        to output additional information. */
    virtual void write();

    /** This function returns true if the instrument records photon packets with the specified
        wavelength, i.e. if the wavelength (redshifted if applicable) falls within the range
        covered by the instrument's wavelength grid. Otherwise it returns false. Although the
        function may return true for a wavelength that falls in a gap between the bins of a grid,
        it never returns false for a wavelength that is recorded. */
    bool recordsWavelength(double lambda) const;

protected:
    /** This function returns the FluxRecorder instance associated with this instrument. This
        function is intended for use in subclasses only. */
    FluxRecorder* instrumentFluxRecorder() { return _recorder; }

    //======================== Data Members =======================

private:
    const WavelengthGrid* _instrumentWavelengthGrid{nullptr};
    FluxRecorder* _recorder{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
