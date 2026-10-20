/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CELLMEDIUM_HPP
#define CELLMEDIUM_HPP

#include "CellMeshInterface.hpp"
#include "ImportedMedium.hpp"

////////////////////////////////////////////////////////////////////

/** A CellMedium instance represents a transfer medium with a spatial density distribution (and,
    optionally, other properties) described by a list of cuboidal cells lined up with the
    coordinate axes; refer to the CellSnapshot class for more information. The cell data is usually
    extracted from a cosmological simulation snapshot, and it must be provided in a column text
    file formatted as described below.

    Refer to the description of the TextInFile class for information on overall formatting and on
    how to include header lines specifying the units for each column in the input file. In case the
    input file has no unit specifications, the default units mentioned below are used instead. The
    number of columns expected in the input file depends on the options configured by the user for
    this CellMedium instance:

    \f[ x_\mathrm{min}\,(\mathrm{pc}) \quad y_\mathrm{min}\,(\mathrm{pc}) \quad
    z_\mathrm{min}\,(\mathrm{pc}) \quad x_\mathrm{max}\,(\mathrm{pc}) \quad
    y_\mathrm{max}\,(\mathrm{pc}) \quad z_\mathrm{max}\,(\mathrm{pc}) \quad \{\,
    \rho\,(\text{M}_\odot\,\text{pc}^{-3}) \;\;|\;\; M\,(\text{M}_\odot) \;\;|\;\;
    n\,(\text{cm}^{-3}) \;\;|\;\; N\,(1) \,\} \f] \f[ [
    v_x\,(\mathrm{km/s}) \quad v_y\,(\mathrm{km/s}) \quad v_z\,(\mathrm{km/s}) ] \quad [
    B_x\,(\mu\mathrm{G}) \quad B_y\,(\mu\mathrm{G}) \quad B_z\,(\mu\mathrm{G}) ] \quad [
    \dots\text{mix params}\dots ] \f]

    The first six columns specify the coordinates of the lower-left and upper-right corners of the
    cell. Depending on the value of the \em massType option, the seventh column lists the average
    mass density \f$\rho\f$, the integrated mass \f$M\f$, the average number density \f$n\f$, or
    the integrated number density \f$N\f$ for the cell.

    If the \em importVelocity option is enabled, the subsequent three columns specify the
    \f$v_x\f$, \f$v_y\f$, \f$v_z\f$ components of the bulk velocity for the material represented
    by the cell.

    If the \em importMagneticField option is enabled, the subsequent three columns specify the
    \f$B_x\f$, \f$B_y\f$, \f$B_z\f$ magnetic field vector components for the cell.

    Finally, the remaining columns specify the quantities requested by the configured material
    mix, such as a temperature, or the parameters used by the configured material mix family to
    select a particular material mix for the cell, in the order listed by the material mix or
    family. */
class CellMedium : public ImportedMedium, public CellMeshInterface
{
    /** The enumeration type indicating the type of mass quantity to be imported. */
    ENUM_DEF(MassType, MassDensity, Mass, NumberDensity, Number)
        ENUM_VAL(MassType, MassDensity, "mass density")
        ENUM_VAL(MassType, Mass, "mass (volume-integrated)")
        ENUM_VAL(MassType, NumberDensity, "number density")
        ENUM_VAL(MassType, Number, "number (volume-integrated)")
    ENUM_END()

    ITEM_CONCRETE(CellMedium, ImportedMedium, "a transfer medium imported from cuboidal cell data")
        ATTRIBUTE_TYPE_INSERT(CellMedium, "CellMeshInterface")

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
        snapshot maintained by this medium. */
    CellSnapshot* cellMesh() const override;

    //===================== Data members ====================

private:
    // an extra pointer to our snapshot used to implement CellMeshInterface (ownership is passed to base class)
    CellSnapshot* _cellSnapshot{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
