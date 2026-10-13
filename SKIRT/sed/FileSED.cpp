/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileSED.hpp"
#include "ColumnInFile.hpp"

//////////////////////////////////////////////////////////////////////

void FileSED::getWavelengthsAndLuminosities(Array& lambdav, Array& pv) const
{
    // read the wavelengths and specific luminosities from the input file
    ColumnInFile infile(this, _filename, "spectral energy distribution");
    infile.addColumn("wavelength", "wavelength", "micron");
    infile.addColumn("specific luminosity", "specific", "W/m");
    infile.readAllColumns(lambdav, pv);
    infile.close();
}

//////////////////////////////////////////////////////////////////////
