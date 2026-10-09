/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CellTreePolicy.hpp"
#include "CellMeshInterface.hpp"
#include "CellSnapshot.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"
#include "TreeNodeEvaluation.hpp"

////////////////////////////////////////////////////////////////////

CellTreePolicy::~CellTreePolicy()
{
    delete _ownMesh;
}

////////////////////////////////////////////////////////////////////

void CellTreePolicy::setupSelfBefore()
{
    TreePolicy::setupSelfBefore();

    if (!filename().empty())
    {
        // import the cells, with a search structure for locating the cell containing a given position
        _ownMesh = new CellSnapshot;
        _ownMesh->open(this, filename(), "cuboidal cells");
        _ownMesh->importBox();
        _ownMesh->setNeedGetEntities();
        _ownMesh->readAndClose();
        _mesh = _ownMesh;
    }
    else
    {
        // search the media and their geometries; don't setup the medium system because we are part of it
        auto ms = find<MediumSystem>(false);
        auto cmi = ms ? ms->interface<CellMeshInterface>(0, 2, false) : nullptr;
        if (!cmi) throw FATALERROR("There is no input file with cells and no medium offering them");
        _mesh = cmi->cellMesh();
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

bool CellTreePolicy::needsSubdivide(TreeNodeEvaluation& node) const
{
    int m = _mesh->cellIndex(Position(node.box().center()));
    return m >= 0 && isCoarser(node.box(), _mesh->boxForCell(m));
}

////////////////////////////////////////////////////////////////////
