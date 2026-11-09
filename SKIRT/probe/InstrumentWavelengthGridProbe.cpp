/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InstrumentWavelengthGridProbe.hpp"
#include "ColumnOutFile.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "Indices.hpp"
#include "Instrument.hpp"
#include "InstrumentSystem.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // returns true if the two wavelength grids have the same bins, so that the files written for them are identical
    bool haveSameBins(const WavelengthGrid* grid1, const WavelengthGrid* grid2)
    {
        if (grid1 == grid2) return true;
        int numBins = grid1->numBins();
        if (grid2->numBins() != numBins) return false;

        // for output units with reverse ordering, the bin width written to the file depends on the type of grid
        bool disjoint1 = dynamic_cast<const DisjointWavelengthGrid*>(grid1) != nullptr;
        bool disjoint2 = dynamic_cast<const DisjointWavelengthGrid*>(grid2) != nullptr;
        if (disjoint1 != disjoint2) return false;

        for (int ell = 0; ell != numBins; ++ell)
        {
            if (grid1->wavelength(ell) != grid2->wavelength(ell) || grid1->leftBorder(ell) != grid2->leftBorder(ell)
                || grid1->rightBorder(ell) != grid2->rightBorder(ell)
                || grid1->effectiveWidth(ell) != grid2->effectiveWidth(ell))
                return false;
        }
        return true;
    }
}

////////////////////////////////////////////////////////////////////

void InstrumentWavelengthGridProbe::probe()
{
    // loop over instruments, skipping an instrument with the same name and the same wavelength grid as an earlier one,
    // such as the SEDInstrument and FrameInstrument that replace a FullInstrument in an upgraded ski file, because
    // the file would be identical; for different grids, the file is written again, so that the conflict is reported
    vector<const Instrument*> written;
    for (auto instrument : find<InstrumentSystem>()->instruments())
    {
        if (std::any_of(written.begin(), written.end(), [instrument](const Instrument* other) {
                return other->instrumentName() == instrument->instrumentName()
                       && haveSameBins(other->instrumentWavelengthGrid(), instrument->instrumentWavelengthGrid());
            }))
            continue;
        written.push_back(instrument);

        writeWavelengthGrid(this, instrument->instrumentWavelengthGrid(), instrument->instrumentName() + "_wavelengths",
                            "wavelengths for instrument " + instrument->instrumentName(),
                            "Wavelength grid for instrument " + instrument->instrumentName());
    }
}

////////////////////////////////////////////////////////////////////

void InstrumentWavelengthGridProbe::writeWavelengthGrid(Probe* item, const WavelengthGrid* wavelengthGrid,
                                                        const string& filename, const string& description,
                                                        const string& longDescription)
{
    auto units = item->find<Units>();

    // create a text file and add the columns
    ColumnOutFile file(item, filename, description);
    file.setLongDescription(longDescription);
    file.addColumn("characteristic wavelength; " + units->swavelength(), units->uwavelength());
    file.addColumn("effective wavelength bin width; " + units->swavelength(), units->uwavelength());
    file.addColumn("left border of wavelength bin; " + units->swavelength(), units->uwavelength());
    file.addColumn("right border of wavelength bin; " + units->swavelength(), units->uwavelength());

    // write the rows
    for (int ell : Indices(wavelengthGrid->numBins(), units->rwavelength()))
    {
        double chara = units->owavelength(wavelengthGrid->wavelength(ell));
        double width = units->owavelength(wavelengthGrid->effectiveWidth(ell));
        double left = units->owavelength(wavelengthGrid->leftBorder(ell));
        double right = units->owavelength(wavelengthGrid->rightBorder(ell));

        // special treatment for output units with reverse ordering
        if (units->rwavelength())
        {
            // swap left and right borders
            std::swap(left, right);

            // fix the bin width
            if (dynamic_cast<const DisjointWavelengthGrid*>(wavelengthGrid))
            {
                // for wavelength grids with constant transmission across the bin,
                // the effective bin width is easily obtained from the borders
                width = right - left;
            }
            else
            {
                // for wavelength grids with variable transmission across the bin,
                // we apply the approximate correction factor dlambda^2 / lambda^2
                double lambda = wavelengthGrid->wavelength(ell);
                double dlambda = wavelengthGrid->effectiveWidth(ell);
                width *= (dlambda * dlambda) / (lambda * lambda);
            }
        }
        file.writeRow(chara, width, left, right);
    }
    file.close();
}

////////////////////////////////////////////////////////////////////
