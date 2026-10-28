/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "IterationHistory.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // returns a description of the specified key for use in error messages
    string keyDescription(const SeriesKey& key)
    {
        string item = key.item ? key.item->type() : "no item";
        return "series " + std::to_string(key.id) + " of " + item;
    }
}

////////////////////////////////////////////////////////////////////

IterationHistory::IterationHistory(SimulationItem* parent)
{
    parent->addChild(this);
}

////////////////////////////////////////////////////////////////////

void IterationHistory::beginLoop(Loop loop)
{
    clear(Lifetime::Loop);
    clear(Lifetime::Grid);
    _loop = loop;
    _iteration = 0;
    _gridChanged = false;
}

////////////////////////////////////////////////////////////////////

void IterationHistory::beginIteration()
{
    if (_gridChanged)
    {
        clear(Lifetime::Grid);
        _gridChanged = false;
    }
    for (auto entry : _scalarOrder) entry->series.advance();
    _iteration++;
}

////////////////////////////////////////////////////////////////////

ScalarSeries& IterationHistory::scalarSeries(const SimulationItem* item, int id, int depth, Lifetime lifetime,
                                             string description, string quantity)
{
    SeriesKey key{item, id};
    auto position = _scalars.find(key);
    if (position != _scalars.end())
    {
        auto& entry = position->second;
        if (entry.lifetime != lifetime)
            throw FATALERROR("The " + keyDescription(key) + " is declared with another lifetime than before");
        if (entry.series.aggregationRule())
            throw FATALERROR("The " + keyDescription(key) + " is declared without the aggregation rule it has");
        entry.series.grow(depth);
        return entry.series;
    }

    auto& entry =
        _scalars.emplace(key, ScalarEntry{ScalarSeries(key, depth, description, quantity), lifetime}).first->second;
    _scalarOrder.push_back(&entry);
    return entry.series;
}

////////////////////////////////////////////////////////////////////

ScalarSeries& IterationHistory::aggregateSeries(const SimulationItem* item, int id, const AggregationRule& rule,
                                                int depth, string description, string quantity)
{
    SeriesKey key{item, id};
    auto position = _scalars.find(key);
    if (position != _scalars.end())
    {
        auto& entry = position->second;
        auto existingRule = entry.series.aggregationRule();
        if (!existingRule || !(*existingRule == rule))
            throw FATALERROR("The " + keyDescription(key) + " is declared with another aggregation rule than before");
        entry.series.grow(depth);
        return entry.series;
    }

    auto& entry =
        _scalars.emplace(key, ScalarEntry{ScalarSeries(key, depth, description, quantity), Lifetime::Simulation})
            .first->second;
    entry.series.setAggregationRule(rule);
    _scalarOrder.push_back(&entry);
    return entry.series;
}

////////////////////////////////////////////////////////////////////

CellWindow& IterationHistory::cellWindow(const SimulationItem* item, int id, int numCells, Lifetime lifetime,
                                         string description)
{
    SeriesKey key{item, id};
    auto position = _windows.find(key);
    if (position != _windows.end())
    {
        auto& entry = position->second;
        if (entry.lifetime != lifetime)
            throw FATALERROR("The cell window " + keyDescription(key)
                             + " is declared with another lifetime than before");
        return entry.window;
    }

    auto& entry = _windows.emplace(key, WindowEntry{CellWindow(key, numCells, description), lifetime}).first->second;
    if (_cellCapacity) entry.window.reserveCells(_cellCapacity);
    _windowOrder.push_back(&entry);
    return entry.window;
}

////////////////////////////////////////////////////////////////////

ScalarSeries& IterationHistory::series(const SimulationItem* item, int id)
{
    SeriesKey key{item, id};
    auto position = _scalars.find(key);
    if (position == _scalars.end()) throw FATALERROR("The " + keyDescription(key) + " has not been declared");
    return position->second.series;
}

////////////////////////////////////////////////////////////////////

vector<ScalarSeries*> IterationHistory::aggregates()
{
    vector<ScalarSeries*> result;
    for (auto entry : _scalarOrder)
        if (entry->series.aggregationRule()) result.push_back(&entry->series);
    return result;
}

////////////////////////////////////////////////////////////////////

void IterationHistory::appendCells(const vector<int>& parentCells)
{
    for (auto entry : _windowOrder) entry->window.appendCells(parentCells);
    _gridChanged = true;
}

////////////////////////////////////////////////////////////////////

void IterationHistory::reserveCells(int numCells)
{
    _cellCapacity = numCells;
    for (auto entry : _windowOrder) entry->window.reserveCells(numCells);
}

////////////////////////////////////////////////////////////////////

vector<const ScalarSeries*> IterationHistory::allScalarSeries() const
{
    vector<const ScalarSeries*> result;
    for (auto entry : _scalarOrder) result.push_back(&entry->series);
    return result;
}

////////////////////////////////////////////////////////////////////

vector<const CellWindow*> IterationHistory::allCellWindows() const
{
    vector<const CellWindow*> result;
    for (auto entry : _windowOrder) result.push_back(&entry->window);
    return result;
}

////////////////////////////////////////////////////////////////////

void IterationHistory::clear(Lifetime lifetime)
{
    for (auto entry : _scalarOrder)
        if (entry->lifetime == lifetime) entry->series.clear();
    for (auto entry : _windowOrder)
        if (entry->lifetime == lifetime) entry->window.reset();
}

////////////////////////////////////////////////////////////////////
