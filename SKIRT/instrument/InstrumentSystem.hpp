/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INSTRUMENTSYSTEM_HPP
#define INSTRUMENTSYSTEM_HPP

#include "Instrument.hpp"

//////////////////////////////////////////////////////////////////////

/** An InstrumentSystem instance keeps a list of zero or more instruments. The instruments can be
    of various nature and do not need to be located at the same observing position. An instrument
    that does not specify its own wavelength grid uses the default wavelength grid of the
    WavelengthGridPool.

    During setup, the instrument system arranges the instruments in groups with the same sight line
    (see Instrument::hasSameSightLine()), regardless of their order in the configuration. The
    simulation launches a single peel-off photon packet for each group, which then serves all
    instruments in the group, so that the peel-off photon packet and the extinction along its path
    are calculated only once. */
class InstrumentSystem : public SimulationItem
{
    ITEM_CONCRETE(InstrumentSystem, SimulationItem, "an instrument system")

        PROPERTY_ITEM_LIST(instruments, Instrument, "the instruments")
        ATTRIBUTE_DEFAULT_VALUE(instruments, "SEDInstrument")
        ATTRIBUTE_REQUIRED_IF(instruments, "false")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function arranges the instruments in groups with the same sight line. Each instrument
        joins the first group whose first instrument has the same sight line, or otherwise starts a
        new group. As a result, the groups are ordered by the position of their first instrument in
        the configuration, and the instruments in each group retain their configuration order. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function flushes any information buffered during photon packet detection for the
        complete instrument system. It calls the flush() function for each of the instruments. */
    void flush();

    /** This function writes the recorded data for the complete instrument system to a set of
        files. It calls the write() function for each of the instruments. */
    void write();

    /** This function returns the groups of instruments with the same sight line, as determined
        during setup. Each group contains at least one instrument, and each instrument in the
        instrument system belongs to exactly one group. */
    const vector<vector<Instrument*>>& sightLineGroups() const { return _sightLineGroups; }

    //======================== Data Members ========================

private:
    // the groups of instruments with the same sight line, initialized during setup
    vector<vector<Instrument*>> _sightLineGroups;
};

////////////////////////////////////////////////////////////////////

#endif
