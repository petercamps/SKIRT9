/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TrivialBoxCriterion.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"
#include "TreeSpatialGrid.hpp"

////////////////////////////////////////////////////////////////////

void TrivialBoxCriterion::setupSelfAfter()
{
    ValueCriterion::setupSelfAfter();

    _box = Box(minX(), minY(), minZ(), maxX(), maxY(), maxZ());
    if (_box.xwidth() <= 0 || _box.ywidth() <= 0 || _box.zwidth() <= 0)
        throw FATALERROR("The box of the refinement criterion should have a positive volume");

    // the medium system and the grid are not set up, because the criterion is part of them
    _ms = find<MediumSystem>(false);
    _grid = find<TreeSpatialGrid>(false);
}

////////////////////////////////////////////////////////////////////

bool TrivialBoxCriterion::changesIn(IterationHistory::Loop loop) const
{
    return loop != IterationHistory::Loop::None;
}

////////////////////////////////////////////////////////////////////

double TrivialBoxCriterion::value(int m) const
{
    return _box.contains(_grid->centralPositionInCell(m)) ? 1. : 0.;
}

////////////////////////////////////////////////////////////////////

double TrivialBoxCriterion::numberDensity(int m) const
{
    double n = 0.;
    int numMedia = _ms->numMedia();
    for (int h = 0; h != numMedia; ++h) n += _ms->numberDensity(m, h);
    return n;
}

////////////////////////////////////////////////////////////////////

string TrivialBoxCriterion::fieldDescription() const
{
    return "position inside the test box";
}

////////////////////////////////////////////////////////////////////
