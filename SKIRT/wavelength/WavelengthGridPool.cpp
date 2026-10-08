/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "WavelengthGridPool.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include <set>

////////////////////////////////////////////////////////////////////

void WavelengthGridPool::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // verify that the names are unique
    std::set<string> names;
    for (auto namedGrid : _wavelengthGrids)
    {
        if (!names.insert(namedGrid->name()).second)
            throw FATALERROR("The wavelength grid name '" + namedGrid->name() + "' occurs more than once in the pool; "
                             + "the pool has " + availableNames());
    }

    // verify that the default grid name, if any, occurs in the list
    defaultGrid();
}

////////////////////////////////////////////////////////////////////

WavelengthGrid* WavelengthGridPool::grid(string name) const
{
    for (auto namedGrid : _wavelengthGrids)
        if (namedGrid->name() == name) return namedGrid->wavelengthGrid();
    throw FATALERROR("There is no wavelength grid named '" + name + "' in the pool; the pool has " + availableNames());
}

////////////////////////////////////////////////////////////////////

WavelengthGrid* WavelengthGridPool::defaultGrid() const
{
    return _defaultGridName.empty() ? nullptr : grid(_defaultGridName);
}

////////////////////////////////////////////////////////////////////

string WavelengthGridPool::availableNames() const
{
    if (_wavelengthGrids.empty()) return "no wavelength grids";

    vector<string> names;
    for (auto namedGrid : _wavelengthGrids) names.push_back("'" + namedGrid->name() + "'");
    return "wavelength grids " + StringUtils::join(names, ", ");
}

////////////////////////////////////////////////////////////////////
