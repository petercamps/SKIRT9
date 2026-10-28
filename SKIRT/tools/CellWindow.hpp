/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CELLWINDOW_HPP
#define CELLWINDOW_HPP

#include "SeriesKey.hpp"

////////////////////////////////////////////////////////////////////

/** A CellWindow holds, for each spatial cell, the running sum of a quantity over the iterations
    accumulated since the last reset, on behalf of a client of the IterationHistory, such as the
    dynamic grid refinement. When cells are subdivided, the history grows the window: a new cell
    starts with the running sum of its parent cell. Instances are created only by the
    IterationHistory class. */
class CellWindow
{
    //============= Construction =============

private:
    // only the IterationHistory class creates and manages cell windows
    friend class IterationHistory;

    /** This constructor creates a window with the specified key, number of cells, and
        description, with all running sums at zero. */
    CellWindow(SeriesKey key, int numCells, string description);

    //============= Accumulating =============

public:
    /** This function adds the value returned by \em valueInCell to the running sum of each cell,
        and increments the number of accumulated iterations. The callback is invoked serially for
        all cells in index order, so that the result is identical on all processes. */
    void accumulate(const std::function<double(int m)>& valueInCell);

    /** This function clears the running sums and the number of accumulated iterations. */
    void reset();

    //============= Querying =============

public:
    /** This function returns the number of iterations accumulated since the last reset. */
    int numIterations() const { return _numIterations; }

    /** This function returns the mean value for the cell with index \em m over the accumulated
        iterations, or NaN if no iterations have been accumulated. */
    double mean(int m) const;

    /** This function returns the number of cells. */
    int numCells() const { return static_cast<int>(_sums.size()); }

    /** This function returns the key of the window. */
    const SeriesKey& key() const { return _key; }

    /** This function returns the human-readable description of the window. */
    const string& description() const { return _description; }

    //============= Managing (by the history) =============

private:
    /** This function grows the window by one cell for each element of \em parentCells. The new
        cells receive consecutive indices following the existing cells, and start with the running
        sum of the indicated parent cell. */
    void appendCells(const vector<int>& parentCells);

    /** This function reserves memory for the specified number of cells, so that the window can grow
        to that number of cells without reallocation. */
    void reserveCells(int numCells);

    //============= Data members =============

private:
    SeriesKey _key;
    string _description;
    vector<double> _sums;   // the running sum for each cell
    int _numIterations{0};  // the number of iterations accumulated since the last reset
};

////////////////////////////////////////////////////////////////////

#endif
