/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileLineSED.hpp"
#include "ColumnInFile.hpp"

//////////////////////////////////////////////////////////////////////

void FileLineSED::getWavelengthsAndLuminosities(Array& lambdav, Array& Lv) const
{
    // read the wavelengths and luminosities from the input file
    ColumnInFile infile(this, _filename, "spectral energy distribution");
    infile.addColumn("wavelength", "wavelength", "micron");
    infile.addColumn("luminosity", "bolluminosity", "Lsun");
    infile.readAllColumns(lambdav, Lv);
    infile.close();
}

//////////////////////////////////////////////////////////////////////
