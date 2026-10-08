/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ReferenceWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "WavelengthGridPool.hpp"

////////////////////////////////////////////////////////////////////

void ReferenceWavelengthGrid::setupSelfBefore()
{
    WavelengthGrid::setupSelfBefore();

    // a pool grid cannot reference another pool grid
    for (Item* ancestor = parent(); ancestor; ancestor = ancestor->parent())
        if (dynamic_cast<WavelengthGridPool*>(ancestor))
            throw FATALERROR("A wavelength grid in the pool cannot reference another pool grid ('" + name() + "')");

    // locate the pool and obtain the referenced grid; setting up the pool also sets up its grids
    auto pool = find<WavelengthGridPool>(false);
    if (!pool)
        throw FATALERROR("There is no wavelength grid pool for the reference to wavelength grid '" + name() + "'");
    pool->setup();
    _grid = pool->grid(name());
}

////////////////////////////////////////////////////////////////////

int ReferenceWavelengthGrid::numBins() const
{
    return _grid->numBins();
}

////////////////////////////////////////////////////////////////////

double ReferenceWavelengthGrid::wavelength(int ell) const
{
    return _grid->wavelength(ell);
}

////////////////////////////////////////////////////////////////////

double ReferenceWavelengthGrid::leftBorder(int ell) const
{
    return _grid->leftBorder(ell);
}

////////////////////////////////////////////////////////////////////

double ReferenceWavelengthGrid::rightBorder(int ell) const
{
    return _grid->rightBorder(ell);
}

////////////////////////////////////////////////////////////////////

double ReferenceWavelengthGrid::effectiveWidth(int ell) const
{
    return _grid->effectiveWidth(ell);
}

////////////////////////////////////////////////////////////////////

double ReferenceWavelengthGrid::transmission(int ell, double lambda) const
{
    return _grid->transmission(ell, lambda);
}

////////////////////////////////////////////////////////////////////

vector<int> ReferenceWavelengthGrid::bins(double lambda) const
{
    return _grid->bins(lambda);
}

////////////////////////////////////////////////////////////////////

int ReferenceWavelengthGrid::bin(double lambda) const
{
    return _grid->bin(lambda);
}

////////////////////////////////////////////////////////////////////

Range ReferenceWavelengthGrid::wavelengthRange() const
{
    return _grid->wavelengthRange();
}

////////////////////////////////////////////////////////////////////
