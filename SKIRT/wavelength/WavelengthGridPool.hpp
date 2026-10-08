/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef WAVELENGTHGRIDPOOL_HPP
#define WAVELENGTHGRIDPOOL_HPP

#include "NamedWavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** A WavelengthGridPool instance holds a list of named wavelength grids that can be used by the
    instruments and probes in the simulation, and designates one of these grids as the default
    wavelength grid. The pool is configured as a property of the MonteCarloSimulation, just before
    the instrument system. It is relevant only for panchromatic simulations, and it is offered
    only to users at the Regular or Expert level; basic users configure a wavelength grid for each
    instrument and probe.

    Each entry in the list is a NamedWavelengthGrid, combining a name with a wavelength grid of any
    type. An instrument or probe can use a pool grid by configuring a ReferenceWavelengthGrid with
    the corresponding name as its wavelength grid, so that a grid can be defined once and used by
    any number of instruments and probes. An instrument or probe without a wavelength grid of its
    own uses the pool grid named by the \em defaultGridName property. If this property is empty,
    there is no default, and each instrument and probe must configure a wavelength grid. The
    default value of this property is empty, because a ski file reader treats an empty string
    attribute as if the attribute were missing, assigning the default value to the property; with a
    nonempty default value, an empty default grid name could not be read back.

    <b>SMILE conditions</b>

    The nonempty list of wavelength grids inserts the global name \c PoolWavelengthGrids, which
    allows configuring a ReferenceWavelengthGrid (see that class for more information). The
    nonempty default grid name inserts the global name \c DefaultInstrumentWavelengthGrid, which
    makes the wavelength grid of instruments and probes optional. The list of wavelength grids
    precedes the default grid name, so that the user first names the grids and then selects the
    default.

    <b>Validation</b>

    Neither the schema nor the configuration wizard can verify the names. During setup, the pool
    therefore reports a fatal error for a duplicate name, and for a default name that does not
    occur in the list. A reference to an unknown name causes a fatal error as well (see the grid()
    function). Each of these error messages lists the available names. */
class WavelengthGridPool : public SimulationItem
{
    ITEM_CONCRETE(WavelengthGridPool, SimulationItem, "a pool of named wavelength grids for instruments and probes")

        PROPERTY_ITEM_LIST(wavelengthGrids, NamedWavelengthGrid, "the named wavelength grids")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthGrids, "NamedWavelengthGrid")
        ATTRIBUTE_INSERT(wavelengthGrids, "wavelengthGrids:PoolWavelengthGrids")

        PROPERTY_STRING(defaultGridName, "the name of the default wavelength grid for instruments and probes")
        ATTRIBUTE_DEFAULT_VALUE(defaultGridName, "")
        ATTRIBUTE_REQUIRED_IF(defaultGridName, "false")
        ATTRIBUTE_INSERT(defaultGridName, "defaultGridName:DefaultInstrumentWavelengthGrid")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the names of the wavelength grids in the pool are unique, and
        that the default grid name, if nonempty, occurs in the list. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the wavelength grid in the pool with the specified name. If there is
        no such grid, the function throws a fatal error listing the available names. The function
        does not require the pool to have been set up, and it does not set up the returned grid. */
    WavelengthGrid* grid(string name) const;

    /** This function returns the default wavelength grid for instruments and probes, i.e. the pool
        grid named by the \em defaultGridName property, or the null pointer if that property is
        empty. If there is no grid with the specified name, the function throws a fatal error
        listing the available names. The function does not require the pool to have been set up,
        and it does not set up the returned grid. */
    WavelengthGrid* defaultGrid() const;

private:
    /** This function returns a list of the names of the wavelength grids in the pool, formatted
        for inclusion in an error message. */
    string availableNames() const;
};

//////////////////////////////////////////////////////////////////////

#endif
