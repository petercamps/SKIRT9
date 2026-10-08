/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef REFERENCEWAVELENGTHGRID_HPP
#define REFERENCEWAVELENGTHGRID_HPP

#include "WavelengthGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** A ReferenceWavelengthGrid instance refers to a named wavelength grid in the WavelengthGridPool
    of the simulation, and forwards all requests to that grid. It can be configured wherever an
    instrument or probe accepts a wavelength grid, so that a grid defined once in the pool can be
    used by any number of instruments and probes. The default grid of the pool can be referenced by
    its name like any other pool grid.

    The instruments and probes do not use the reference itself. Instead, they obtain the
    referenced grid through the Configuration::wavelengthGrid() function, so that they can take
    advantage of the specific type of that grid (e.g. for a grid with disjoint bins). The
    forwarding functions serve the other clients, such as the functions that determine the
    wavelength range of the simulation.

    <b>SMILE conditions</b>

    This type is allowed only if the pool has at least one wavelength grid, which inserts the name
    \c PoolWavelengthGrids, and if the instrument system or the probe system has been added to the
    configuration, which inserts the name \c InstrumentSystem or \c ProbeSystem. Each instrument and
    probe is configured inside one of these systems, which follow the pool, so that the second
    condition ensures that the type is offered for instruments and probes, but not for a grid
    inside the pool, which excludes reference cycles. Both system names are needed, because a ski
    file may omit the instrument system, which then is added with its default value only after
    the probe system has been read. The first condition by itself does not suffice, because the
    name inserted by a nonempty item list is already in effect while the items in that list are
    being configured.

    <b>Validation</b>

    During setup, a reference reports a fatal error if the pool does not have a grid with the
    configured name, listing the available names. It also reports a fatal error if it is located
    inside the pool, which can happen only in a manually edited ski file. */
class ReferenceWavelengthGrid : public WavelengthGrid
{
    ITEM_CONCRETE(ReferenceWavelengthGrid, WavelengthGrid, "a reference to a named wavelength grid in the pool")
        ATTRIBUTE_TYPE_ALLOWED_IF(ReferenceWavelengthGrid, "PoolWavelengthGrids&(InstrumentSystem|ProbeSystem)")

        PROPERTY_STRING(name, "the name of the referenced wavelength grid in the pool")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function locates the wavelength grid pool, sets it up, and obtains the referenced
        grid. It throws a fatal error if the reference is located inside the pool, if there is no
        pool, or if the pool does not have a grid with the configured name. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the referenced wavelength grid in the pool. It may be called only
        after setup has been completed. */
    WavelengthGrid* referencedGrid() const { return _grid; }

    //======================== Forwarding =======================

    /** This function returns the number of bins in the referenced grid. */
    int numBins() const override;

    /** This function returns the characteristic wavelength for the specified bin in the
        referenced grid. */
    double wavelength(int ell) const override;

    /** This function returns the left border of the specified bin in the referenced grid. */
    double leftBorder(int ell) const override;

    /** This function returns the right border of the specified bin in the referenced grid. */
    double rightBorder(int ell) const override;

    /** This function returns the effective width of the specified bin in the referenced grid. */
    double effectiveWidth(int ell) const override;

    /** This function returns the relative transmission for the specified bin in the referenced
        grid at the specified wavelength. */
    double transmission(int ell, double lambda) const override;

    /** This function returns the indices of the bins in the referenced grid that may have a
        nonzero transmission at the specified wavelength. */
    vector<int> bins(double lambda) const override;

    /** This function returns the index of one of the bins in the referenced grid that may have a
        nonzero transmission at the specified wavelength, or -1 if there is no such bin. */
    int bin(double lambda) const override;

    /** This function returns the wavelength range covered by the referenced grid. */
    Range wavelengthRange() const override;

    //======================== Data Members ========================

private:
    WavelengthGrid* _grid{nullptr};  // the referenced grid, initialized during setup
};

//////////////////////////////////////////////////////////////////////

#endif
