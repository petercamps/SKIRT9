/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TreeSpatialGridTopologyProbe.hpp"
#include "ColumnOutFile.hpp"
#include "MediumSystem.hpp"
#include "TreeSpatialGrid.hpp"

////////////////////////////////////////////////////////////////////

void TreeSpatialGridTopologyProbe::probe()
{
    // locate the grid (it is OK for the medium system to have no media components)
    auto ms = find<MediumSystem>(false);
    if (auto grid = ms ? ms->find<TreeSpatialGrid>(false) : nullptr)
    {
        ColumnOutFile outfile(this, itemName() + "_treetop", "spatial tree grid topology");
        grid->writeTopology(&outfile);
        outfile.close();
    }
}

////////////////////////////////////////////////////////////////////
