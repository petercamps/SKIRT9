/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TreeSpatialGridTopologyProbe.hpp"
#include "ColumnOutFile.hpp"
#include "MediumSystem.hpp"
#include "TreeSpatialGrid.hpp"

////////////////////////////////////////////////////////////////////

Probe::When TreeSpatialGridTopologyProbe::when() const
{
    switch (probeAfter())
    {
        case ProbeAfter::Setup: return When::Setup;
        case ProbeAfter::Run: return When::Run;
        case ProbeAfter::Primary: return When::Primary;
        case ProbeAfter::Secondary: return When::Secondary;
    }
    return When::Setup;
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGridTopologyProbe::probe()
{
    // locate the grid (it is OK for the medium system to have no media components)
    auto ms = find<MediumSystem>(false);
    if (auto grid = ms ? ms->find<TreeSpatialGrid>(false) : nullptr)
    {
        // after an iteration, include the zero-padded iteration index in the file name
        string prefix = itemName();
        if (iter() > 0)
        {
            string index = std::to_string(iter());
            if (index.length() < 3) index.insert(0, 3 - index.length(), '0');
            prefix += index;
        }

        ColumnOutFile outfile(this, prefix + "_treetop", "spatial tree grid topology");
        grid->writeTopology(&outfile);
        outfile.close();
    }
}

////////////////////////////////////////////////////////////////////
