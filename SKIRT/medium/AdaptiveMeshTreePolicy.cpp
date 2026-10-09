/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshTreePolicy.hpp"
#include "AdaptiveMeshInterface.hpp"
#include "AdaptiveMeshSnapshot.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"
#include "TreeNodeEvaluation.hpp"
#include "TreeSpatialGrid.hpp"

////////////////////////////////////////////////////////////////////

AdaptiveMeshTreePolicy::~AdaptiveMeshTreePolicy()
{
    delete _ownMesh;
}

////////////////////////////////////////////////////////////////////

void AdaptiveMeshTreePolicy::setupSelfBefore()
{
    TreePolicy::setupSelfBefore();

    if (!filename().empty())
    {
        // import the mesh, using the domain of the grid; don't setup the grid because we are part of it
        _ownMesh = new AdaptiveMeshSnapshot;
        _ownMesh->open(this, filename(), "adaptive mesh");
        _ownMesh->setExtent(find<TreeSpatialGrid>(false)->extent());
        _ownMesh->readAndClose();
        _mesh = _ownMesh;
    }
    else
    {
        // search the media and their geometries; don't setup the medium system because we are part of it
        auto ms = find<MediumSystem>(false);
        auto ami = ms ? ms->interface<AdaptiveMeshInterface>(0, 2, false) : nullptr;
        if (!ami) throw FATALERROR("There is no input file with an adaptive mesh and no medium offering one");
        _mesh = ami->adaptiveMesh();
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns true if the node is wider than the cell along any axis, allowing for rounding errors
    bool isCoarser(const Box& node, const Box& cell)
    {
        const double factor = 1. + 1e-9;
        return node.xwidth() > factor * cell.xwidth() || node.ywidth() > factor * cell.ywidth()
               || node.zwidth() > factor * cell.zwidth();
    }
}

////////////////////////////////////////////////////////////////////

bool AdaptiveMeshTreePolicy::needsSubdivide(TreeNodeEvaluation& node) const
{
    int m = _mesh->cellIndex(Position(node.box().center()));
    return m >= 0 && isCoarser(node.box(), _mesh->extent(m));
}

////////////////////////////////////////////////////////////////////
