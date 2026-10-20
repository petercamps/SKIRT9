/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CELLGEOMETRY_HPP
#define CELLGEOMETRY_HPP

#include "CellMeshInterface.hpp"
#include "ImportedGeometry.hpp"

////////////////////////////////////////////////////////////////////

/** A CellGeometry instance represents a 3D geometry with a spatial density distribution described
    by a list of cuboidal cells lined up with the coordinate axes; refer to the CellSnapshot class
    for more information. The cell data is usually extracted from a cosmological simulation
    snapshot, and it must be provided in a column text file formatted as described below. The total
    mass in the geometry is normalized to unity after importing the data.

    Refer to the description of the TextInFile class for information on overall formatting and on
    how to include header lines specifying the units for each column in the input file. In case the
    input file has no unit specifications, the default units mentioned below are used instead. The
    input file should contain 7 columns:

    \f[ x_\mathrm{min}\,(\mathrm{pc}) \quad y_\mathrm{min}\,(\mathrm{pc}) \quad
    z_\mathrm{min}\,(\mathrm{pc}) \quad x_\mathrm{max}\,(\mathrm{pc}) \quad
    y_\mathrm{max}\,(\mathrm{pc}) \quad z_\mathrm{max}\,(\mathrm{pc}) \quad \{\,
    \rho\,(\text{M}_\odot\,\text{pc}^{-3}) \;\;|\;\; M\,(\text{M}_\odot) \;\;|\;\;
    n\,(\text{cm}^{-3}) \;\;|\;\; N\,(1) \,\} \f]

    The first six columns specify the coordinates of the lower-left and upper-right corners of the
    cell. Depending on the value of the \em massType option, the seventh column lists the average
    mass density \f$\rho\f$, the integrated mass \f$M\f$, the average number density \f$n\f$, or
    the integrated number density \f$N\f$ for the cell. The precise units for this field are
    irrelevant because the total mass in the geometry will be normalized to unity after importing
    the data. However, the import procedure still insists on knowing the units. */
class CellGeometry : public ImportedGeometry, public CellMeshInterface
{
    /** The enumeration type indicating the type of mass quantity to be imported. */
    ENUM_DEF(MassType, MassDensity, Mass, NumberDensity, Number)
        ENUM_VAL(MassType, MassDensity, "mass density")
        ENUM_VAL(MassType, Mass, "mass (volume-integrated mass density)")
        ENUM_VAL(MassType, NumberDensity, "number density")
        ENUM_VAL(MassType, Number, "number (volume-integrated number density)")
    ENUM_END()

    ITEM_CONCRETE(CellGeometry, ImportedGeometry, "a geometry imported from cuboidal cell data")
        ATTRIBUTE_TYPE_INSERT(CellGeometry, "CellMeshInterface")

        PROPERTY_ENUM(massType, MassType, "the type of mass quantity to be imported")
        ATTRIBUTE_DEFAULT_VALUE(massType, "MassDensity")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new CellSnapshot object, calls its open() function, configures
        it to import a mass or density column, and finally returns a pointer to the object.
        Ownership of the Snapshot object is transferred to the caller. */
    Snapshot* createAndOpenSnapshot() override;

    //=================== Other functions ==================

protected:
    /** This function implements the CellMeshInterface interface. It returns a pointer to the cell
        snapshot maintained by this geometry. */
    CellSnapshot* cellMesh() const override;

    //===================== Data members ====================

private:
    // an extra pointer to our snapshot used to implement CellMeshInterface (ownership is passed to base class)
    CellSnapshot* _cellSnapshot{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
