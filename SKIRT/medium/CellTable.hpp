/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CELLTABLE_HPP
#define CELLTABLE_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** A CellTable instance holds a two-dimensional table of double values with a row for each spatial
    cell and a fixed number of columns, for example a radiation field table with a column for each
    wavelength bin. The values are stored row by row in a standard vector, so that capacity for
    additional rows can be reserved in advance, and the table can grow with the spatial grid during
    dynamic grid refinement without reallocation. */
class CellTable
{
public:
    /** This function resizes the table to the specified number of rows (cells) and columns, and
        sets all values to zero. */
    void resize(int numCells, int numColumns)
    {
        _numColumns = numColumns;
        _data.assign(static_cast<size_t>(numCells) * _numColumns, 0.);
    }

    /** This function reserves memory for the specified number of rows (cells), so that the table
        can grow to that number of rows without reallocation. */
    void reserveCells(int numCells) { _data.reserve(static_cast<size_t>(numCells) * _numColumns); }

    /** This function grows the table by one row for each element of \em parentCells. The new rows
        receive consecutive indices following the existing rows, and start with the values of the
        indicated parent row. */
    void appendCells(const vector<int>& parentCells)
    {
        size_t numOld = numCells();
        _data.resize((numOld + parentCells.size()) * _numColumns);
        for (size_t i = 0; i != parentCells.size(); ++i)
        {
            auto parent = _data.begin() + parentCells[i] * _numColumns;
            std::copy(parent, parent + _numColumns, _data.begin() + (numOld + i) * _numColumns);
        }
    }

    /** This function multiplies the values in the row with index \f$m\f$ by the specified factor.
        */
    void scaleCell(int m, double factor)
    {
        auto row = _data.begin() + m * _numColumns;
        for (auto value = row; value != row + _numColumns; ++value) *value *= factor;
    }

    /** This function sets all values in the table to zero. */
    void setToZero() { std::fill(_data.begin(), _data.end(), 0.); }

    /** This function returns the number of rows (cells) in the table. */
    int numCells() const { return _numColumns ? static_cast<int>(_data.size() / _numColumns) : 0; }

    /** This function returns the total number of values in the table. */
    size_t size() const { return _data.size(); }

    /** This function returns the number of bytes held by the table for each row (cell). */
    size_t bytesPerCell() const { return _numColumns * sizeof(double); }

    /** This function returns the value in the specified row (cell) and column. */
    double operator()(int m, int ell) const { return _data[m * _numColumns + ell]; }

    /** This function returns a writable reference to the value in the specified row (cell) and
        column. */
    double& operator()(int m, int ell) { return _data[m * _numColumns + ell]; }

    /** This function returns a writable reference to the vector holding the values, row by row,
        for example for communicating the values between processes. The size of the vector must
        not be changed through this reference. */
    vector<double>& data() { return _data; }

private:
    vector<double> _data;
    size_t _numColumns{0};
};

////////////////////////////////////////////////////////////////////

#endif
