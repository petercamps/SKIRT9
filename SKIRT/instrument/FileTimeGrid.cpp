/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileTimeGrid.hpp"
#include "ColumnInFile.hpp"

////////////////////////////////////////////////////////////////////

void FileTimeGrid::getTimeBins(vector<Bin>& bins) const
{
    // setup the columns for the input file
    ColumnInFile infile(this, _filename, "time grid");
    infile.addColumn("time", "timelag", "s");
    infile.addColumn("left", "timelag", "s");
    infile.addColumn("right", "timelag", "s");

    // read the time bins from the input file
    while (true)
    {
        double time, left, right;
        if (!infile.readRow(time, left, right)) break;
        bins.emplace_back(left, time, right);
    }
    infile.close();
}

//////////////////////////////////////////////////////////////////////
