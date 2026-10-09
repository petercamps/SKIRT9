/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CELLMESHINTERFACE_HPP
#define CELLMESHINTERFACE_HPP

#include "Basics.hpp"
class CellSnapshot;

////////////////////////////////////////////////////////////////////

/** CellMeshInterface is a pure interface. It provides access to the cell snapshot (a list of
    cuboidal cells lined up with the coordinate axes) maintained by the object that implements the
    interface. */
class CellMeshInterface
{
protected:
    /** The empty constructor for the interface. */
    CellMeshInterface() {}

public:
    /** The empty destructor for the interface. */
    virtual ~CellMeshInterface() {}

    /** This function must be implemented in a derived class. It returns a pointer to the cell
        snapshot maintained by the object that implements the interface. */
    virtual CellSnapshot* cellMesh() const = 0;
};

/////////////////////////////////////////////////////////////////////////////

#endif
