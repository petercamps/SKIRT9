/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileMesh.hpp"
#include "ColumnInFile.hpp"

////////////////////////////////////////////////////////////////////

vector<double> FileMesh::getMeshBorderPoints() const
{
    // read the mesh border points from the input file
    vector<double> points;
    ColumnInFile infile(this, _filename, "mesh border points");
    infile.addColumn("border point");
    Array row;
    while (infile.readRow(row)) points.push_back(row[0]);
    infile.close();

    return points;
}

////////////////////////////////////////////////////////////////////
