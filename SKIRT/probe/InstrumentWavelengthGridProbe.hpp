/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WAVELENGTHGRIDPROBE_HPP
#define WAVELENGTHGRIDPROBE_HPP

#include "SpecialtyProbe.hpp"
class WavelengthGrid;

////////////////////////////////////////////////////////////////////

/** InstrumentWavelengthGridProbe outputs a text column file with information on the wavelength grid
    for each instrument. Each file is named <tt>prefix_instr_wavelengths.dat</tt> so that it sits
    next to the files written by the corresponding instrument. For each wavelength bin, the file
    lists the characteristic wavelength, the wavelength bin width, and the left and right borders of
    the bin.

    Several instruments may have the same name, as long as they write output files with different
    names, as do an SEDInstrument and a FrameInstrument. If such instruments have the same
    wavelength grid, the probe writes a single file for them. If their grids differ, it writes the
    file for each of them, so that only the last one is kept, and a warning about the conflicting
    file names is logged (see the FilePaths class). */
class InstrumentWavelengthGridProbe : public SpecialtyProbe
{
    ITEM_CONCRETE(InstrumentWavelengthGridProbe, SpecialtyProbe, "wavelength grid: instruments")
        ATTRIBUTE_TYPE_DISPLAYED_IF(InstrumentWavelengthGridProbe, "Level2&Instrument")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;

public:
    /** This function outputs a column text file for the specified wavelength grid in the format as
        described in the header of this class. It can be used from other probes to output
        wavelength grid details. The short \em description is used in the log message, and the
        \em longDescription is written as the first header line of the file. */
    static void writeWavelengthGrid(Probe* item, const WavelengthGrid* wavelengthGrid, const string& filename,
                                    const string& description, const string& longDescription);
};

////////////////////////////////////////////////////////////////////

#endif
