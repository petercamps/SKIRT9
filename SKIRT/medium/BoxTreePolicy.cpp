/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BoxTreePolicy.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "TreeNodeEvaluation.hpp"
#include "TreeSpatialGrid.hpp"

////////////////////////////////////////////////////////////////////

void BoxTreePolicy::setupSelfBefore()
{
    TreePolicy::setupSelfBefore();

    // verify that the box is not empty
    if (_maxX <= _minX) throw FATALERROR("The extent of the box should be positive in the X direction");
    if (_maxY <= _minY) throw FATALERROR("The extent of the box should be positive in the Y direction");
    if (_maxZ <= _minZ) throw FATALERROR("The extent of the box should be positive in the Z direction");

    // copy the coordinates to a Box instance for ease of use
    _box = Box(_minX, _minY, _minZ, _maxX, _maxY, _maxZ);

    // verify that the box is inside the spatial grid; don't setup the grid because we are part of it
    auto grid = find<TreeSpatialGrid>(false);
    if (!grid->extent().contains(_box))
        find<Log>()->warning("The box of the " + type() + " is not fully inside the spatial grid domain");
}

////////////////////////////////////////////////////////////////////

bool BoxTreePolicy::needsSubdivide(TreeNodeEvaluation& node) const
{
    return _box.intersects(node.box()) && _policy->needsSubdivide(node);
}

////////////////////////////////////////////////////////////////////
