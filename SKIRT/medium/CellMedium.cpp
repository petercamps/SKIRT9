/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CellMedium.hpp"
#include "CellSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* CellMedium::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new CellSnapshot;
    _cellSnapshot = snapshot;
    snapshot->open(this, filename(), "cuboidal cells");

    // honor custom column reordering
    snapshot->useColumns(useColumns());

    // configure the cell bounding box columns
    snapshot->importBox();

    // configure the mass or density column
    switch (massType())
    {
        case MassType::MassDensity: snapshot->importMassDensity(); break;
        case MassType::Mass: snapshot->importMass(); break;
        case MassType::NumberDensity: snapshot->importNumberDensity(); break;
        case MassType::Number: snapshot->importNumber(); break;
    }
    return snapshot;
}

////////////////////////////////////////////////////////////////////

CellSnapshot* CellMedium::cellMesh() const
{
    return _cellSnapshot;
}

////////////////////////////////////////////////////////////////////
