/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PEELOFFINSTRUMENT_HPP
#define PEELOFFINSTRUMENT_HPP

#include "Direction.hpp"
#include "Instrument.hpp"
#include "Position.hpp"

////////////////////////////////////////////////////////////////////

/** PeelOffInstrument is an abstract class representing instruments that record peel-off photon
    packets. At each emission and scattering event, the simulation sends a peel-off photon packet
    towards each distinct sight line, and the instruments on that sight line record it after
    attenuation along the path. Subclasses implement instruments that can vary in type of
    projection (e.g. parallel or perspective) and in what is recorded (e.g. an SED or a data cube).
    This class offers the interface through which the simulation determines where a peel-off photon
    packet is detected, and the user property for recording polarization, which is meaningful only
    for peel-off photon packets.

    The simulation detects a peel-off photon packet in a number of steps, so that it can avoid
    unnecessary work. Firstly, the instrument system groups instruments with the same sight line,
    using the hasSameSightLine() function, so that a single peel-off photon packet can serve all
    instruments in a group. For each instrument, the simulation then verifies that the instrument
    records the peel-off photon packet's wavelength, using the recordsWavelength() function, and
    obtains the pixel hit by the photon packet (if any), using the locate() function. Only if the
    photon packet is actually detected by at least one instrument in the group, the simulation
    calculates the extinction along the sight line and passes it to the detect() function of the
    instruments concerned. */
class PeelOffInstrument : public Instrument
{
    ITEM_ABSTRACT(PeelOffInstrument, Instrument, "an instrument that records peel-off photon packets")

        ATTRIBUTE_SUB_PROPERTIES_HERE(PeelOffInstrument)

        PROPERTY_BOOL(recordPolarization, "record polarization (Stokes vector elements)")
        ATTRIBUTE_DEFAULT_VALUE(recordPolarization, "false")
        ATTRIBUTE_DISPLAYED_IF(recordPolarization, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function performs setup for the instrument. In addition to the setup performed by the
        base class, it passes the value of the \em recordPolarization property to the flux
        recorder. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This structure describes where an instrument detects a photon packet, as determined by the
        instrument's geometry and returned by the locate() function. */
    struct Detection
    {
        /** The index of the pixel in the instrument frame where the photon packet is detected,
            zero for instruments without spatial resolution, or -1 if the photon packet is not
            detected. */
        int pixel{-1};

        /** The distance from the photon packet's launching position to the instrument, used to
            calibrate the flux and to limit the path along which the extinction is calculated.
            This distance is infinite for distant instruments. */
        double distance{std::numeric_limits<double>::infinity()};
    };

    /** This function records the detection of the specified photon packet by the instrument, at
        the pixel and distance determined by the locate() function, and with the specified
        extinction factor \f$\exp(-\tau)\f$ (i.e. 1 for a simulation without media). It simply
        calls the corresponding function of the FluxRecorder instance associated with this
        instrument. */
    void detect(const PhotonPacket* pp, Detection detection, double extinction);

    //=========== Functions to be implemented in subclass ===========

public:
    /** This function returns true if the specified instrument has the same sight line as the
        receiving instrument. Two instruments have the same sight line if a peel-off photon packet
        launched towards one of them, including the extinction calculated along its path, can be
        used for the other one as well. In other words, for any launching position, both
        instruments must have the same direction towards the observer (see bfkobs()) and, for local
        instruments, the same distance to the observer. If any of the two instruments records
        polarization, they must also have the same orientation of the instrument frame (see
        bfky()). The InstrumentSystem calls this function during setup to group instruments with
        the same sight line. The implementation must be provided in a subclass. */
    virtual bool hasSameSightLine(const PeelOffInstrument* other) const = 0;

    /** This function returns the direction towards the observer, expressed in model coordinates,
        given the photon packet's launching position. The implementation must be provided in a
        subclass. */
    virtual Direction bfkobs(Position bfr) const = 0;

    /** This function returns the direction along the positive y-axis of the instrument frame,
        expressed in model coordinates, given the photon packet's launching position. The
        implementation must be provided in a subclass. */
    virtual Direction bfky(Position bfr) const = 0;

    /** This function determines where the instrument would detect a peel-off photon packet
        launched towards it from the specified position. It returns the pixel index and distance
        as described for the Detection structure, with a pixel index of -1 if the photon packet
        would not be detected, e.g. because it projects outside of the instrument's aperture or
        field of view. The implementation must be provided in a subclass. */
    virtual Detection locate(Position bfr) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
