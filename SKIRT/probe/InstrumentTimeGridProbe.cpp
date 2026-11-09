/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InstrumentTimeGridProbe.hpp"
#include "ColumnOutFile.hpp"
#include "InstrumentSystem.hpp"
#include "TimeGrid.hpp"
#include "TimeInstrument.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // returns true if the two time grids have the same bins, so that the files written for them are identical
    bool haveSameBins(const TimeGrid* grid1, const TimeGrid* grid2)
    {
        if (grid1 == grid2) return true;
        int numBins = grid1->numBins();
        if (grid2->numBins() != numBins) return false;
        for (int k = 0; k != numBins; ++k)
        {
            if (grid1->time(k) != grid2->time(k) || grid1->left(k) != grid2->left(k)
                || grid1->right(k) != grid2->right(k) || grid1->width(k) != grid2->width(k))
                return false;
        }
        return true;
    }
}

////////////////////////////////////////////////////////////////////

void InstrumentTimeGridProbe::probe()
{
    // loop over time instruments, skipping an instrument with the same name and the same time grid as an earlier one,
    // because the file would be identical; for different grids, the file is written again, so that the conflict is
    // reported
    vector<const TimeInstrument*> written;
    for (auto instrument : find<InstrumentSystem>()->instruments())
    {
        if (auto timeInstrument = dynamic_cast<TimeInstrument*>(instrument))
        {
            if (std::any_of(written.begin(), written.end(), [timeInstrument](const TimeInstrument* other) {
                    return other->instrumentName() == timeInstrument->instrumentName()
                           && haveSameBins(other->timeGrid(), timeInstrument->timeGrid());
                }))
                continue;
            written.push_back(timeInstrument);

            auto units = find<Units>();

            // create a text file and add the columns
            string filename = timeInstrument->instrumentName() + "_timegrid";
            string description = "time grid for instrument " + timeInstrument->instrumentName();
            ColumnOutFile file(this, filename, description);
            file.setLongDescription("Time grid for instrument " + timeInstrument->instrumentName());
            file.addColumn("characteristic time", units->utimelag());
            file.addColumn("left border of time bin", units->utimelag());
            file.addColumn("right border of time bin", units->utimelag());
            file.addColumn("width of time bin", units->utimelag());

            // write the rows
            auto timeGrid = timeInstrument->timeGrid();
            int numBins = timeGrid->numBins();
            for (int k = 0; k != numBins; ++k)
            {
                double time = units->otimelag(timeGrid->time(k));
                double left = units->otimelag(timeGrid->left(k));
                double right = units->otimelag(timeGrid->right(k));
                double width = units->otimelag(timeGrid->width(k));
                file.writeRow(time, left, right, width);
            }
            file.close();
        }
    }
}

////////////////////////////////////////////////////////////////////
