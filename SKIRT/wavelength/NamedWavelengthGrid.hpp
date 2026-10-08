/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef NAMEDWAVELENGTHGRID_HPP
#define NAMEDWAVELENGTHGRID_HPP

#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** A NamedWavelengthGrid instance combines a name with a wavelength grid of any type. It serves as
    an entry in the list of wavelength grids held by the WavelengthGridPool, so that wavelength
    grids outside of the pool do not need to carry a name. The name identifies the wavelength grid
    within the pool; it is used to designate the default wavelength grid for instruments and probes
    and to reference the wavelength grid through a ReferenceWavelengthGrid. */
class NamedWavelengthGrid : public SimulationItem
{
    ITEM_CONCRETE(NamedWavelengthGrid, SimulationItem, "a named wavelength grid")

        PROPERTY_STRING(name, "the name of the wavelength grid")
        ATTRIBUTE_DEFAULT_VALUE(name, "default")

        PROPERTY_ITEM(wavelengthGrid, WavelengthGrid, "the wavelength grid")
        ATTRIBUTE_DEFAULT_VALUE(wavelengthGrid, "LogWavelengthGrid")

    ITEM_END()
};

//////////////////////////////////////////////////////////////////////

#endif
