/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ATPOSITIONSFORM_HPP
#define ATPOSITIONSFORM_HPP

#include "Array.hpp"
#include "GenericForm.hpp"

//////////////////////////////////////////////////////////////////////

/** AtPositionsForm represents a generic probe form. Refer to the ProbeFormBridge class for more
    information about probes and forms.

    This particular form outputs a text column file listing the quantity being probed at the
    positions listed in a user-provided input file. This allows, for example, probing at the
    positions in a Particle or VoronoiMesh import file used for defining sources or media by
    configuring this probe with the same import file (any columns other than the position
    coordinates will be ignored).

    The input text column file is expected to have three columns, respectively representing the x,
    y and z coordinates of a position in the spatial domain of the simulation. In case the input
    file has no unit specifications, the default units are parsec. Refer to the description of the
    TextInFile class for more information on overall formatting, on how to include header lines
    specifying the units for each column, and on how to use the \em useColumns property to address
    columns other than the first three.

    The output file contains a line for each (non-empty and non-comment) line in the input file, in
    the same order. The first three columns repeat the x, y and z coordinates of the position (in
    the configured output units, which may differ from the input units). Subsequent column(s) list
    the quantity being probed.

    The positions are read from the input file during setup and kept in memory, because a
    simulation reads all of its input files during setup (see the I/O cutover described for the
    FilePaths class). The input file can also be provided as a bundle in an HDF5 file (see the
    ColumnInFile class). */
class AtPositionsForm : public GenericForm
{
    ITEM_CONCRETE(AtPositionsForm, GenericForm, "a text column file with values for each imported position")
        ATTRIBUTE_TYPE_DISPLAYED_IF(AtPositionsForm, "Level2")

        PROPERTY_STRING(filename, "the name of the file listing the positions")

        PROPERTY_STRING(useColumns, "a list of names corresponding to columns in the file to be imported")
        ATTRIBUTE_DEFAULT_VALUE(useColumns, "")
        ATTRIBUTE_REQUIRED_IF(useColumns, "false")
        ATTRIBUTE_DISPLAYED_IF(useColumns, "Level3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads the positions from the input file. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function causes the form to output file(s) as described in the class header for the
        quantity being probed according to the information provided by the specified
        ProbeFormBridge instance. */
    void writeQuantity(const ProbeFormBridge* bridge) const override;

    //======================== Data Members ========================

private:
    // the coordinates of the positions read from the input file, in internal units
    Array _xv;
    Array _yv;
    Array _zv;
};

//////////////////////////////////////////////////////////////////////

#endif
