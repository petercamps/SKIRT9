/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INSTRUMENTTIMEGRIDPROBE_HPP
#define INSTRUMENTTIMEGRIDPROBE_HPP

#include "SpecialtyProbe.hpp"
class WavelengthGrid;

////////////////////////////////////////////////////////////////////

/** InstrumentTimeGridProbe outputs a text column file with information on the time grid for each
    time instrument in the simulation. Each file is named <tt>prefix_instr_timegrid.dat</tt> so
    that it sits next to the files written by the corresponding instrument. For each time bin, the
    file lists the characteristic time, the left border, the right border, and the width of the
    bin.

    If several time instruments with the same name have the same time grid, the probe writes a
    single file for them. If their grids differ, it writes the file for each of them, so that only
    the last one is kept, and a warning about the conflicting file names is logged (see the
    FilePaths class). */
class InstrumentTimeGridProbe : public SpecialtyProbe
{
    ITEM_CONCRETE(InstrumentTimeGridProbe, SpecialtyProbe, "time grid: instruments")
        ATTRIBUTE_TYPE_DISPLAYED_IF(InstrumentTimeGridProbe, "Level2&TimeInstrument")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
