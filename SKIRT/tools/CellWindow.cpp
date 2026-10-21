/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CellWindow.hpp"

////////////////////////////////////////////////////////////////////

CellWindow::CellWindow(SeriesKey key, int numCells, string description)
    : _key{key}, _description{description}, _sums(max(numCells, 0))
{}

////////////////////////////////////////////////////////////////////

void CellWindow::accumulate(const std::function<double(int m)>& valueInCell)
{
    int n = numCells();
    for (int m = 0; m != n; ++m) _sums[m] += valueInCell(m);
    _numIterations++;
}

////////////////////////////////////////////////////////////////////

void CellWindow::reset()
{
    std::fill(_sums.begin(), _sums.end(), 0.);
    _numIterations = 0;
}

////////////////////////////////////////////////////////////////////

double CellWindow::mean(int m) const
{
    return _numIterations ? _sums[m] / _numIterations : std::numeric_limits<double>::quiet_NaN();
}

////////////////////////////////////////////////////////////////////

void CellWindow::appendCells(const vector<int>& parentCells)
{
    _sums.reserve(_sums.size() + parentCells.size());
    for (int parent : parentCells) _sums.push_back(_sums[parent]);
}

////////////////////////////////////////////////////////////////////
