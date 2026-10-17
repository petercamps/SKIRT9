/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AtPositionsForm.hpp"
#include "ColumnInFile.hpp"
#include "ColumnOutFile.hpp"
#include "ProbeFormBridge.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

void AtPositionsForm::setupSelfBefore()
{
    GenericForm::setupSelfBefore();

    // read the positions from the input file
    ColumnInFile infile(this, filename(), "positions");
    infile.useColumns(useColumns());
    infile.addColumn("position x", "length", "pc");
    infile.addColumn("position y", "length", "pc");
    infile.addColumn("position z", "length", "pc");
    infile.readAllColumns(_xv, _yv, _zv);
}

////////////////////////////////////////////////////////////////////

void AtPositionsForm::writeQuantity(const ProbeFormBridge* bridge) const
{
    auto units = bridge->units();

    // create a text column file and add the column definitions
    ColumnOutFile outfile(bridge->probe(), bridge->prefix(), bridge->description() + " at imported positions");
    outfile.setLongDescription(StringUtils::toUpperFirst(bridge->description()) + " at imported positions");
    outfile.addColumn("position x", units->ulength());
    outfile.addColumn("position y", units->ulength());
    outfile.addColumn("position z", units->ulength());
    bridge->addColumnDefinitions(outfile);

    // allocate room for a single row
    int numValues = bridge->numValues();
    Array values(numValues);
    vector<double> columns(numValues + 3);

    // write a line for each position in the input file
    size_t numPositions = _xv.size();
    for (size_t i = 0; i != numPositions; ++i)
    {
        double x = _xv[i];
        double y = _yv[i];
        double z = _zv[i];
        bridge->valuesAtPosition(Position(x, y, z), values);
        columns[0] = units->olength(x);
        columns[1] = units->olength(y);
        columns[2] = units->olength(z);
        for (int p = 0; p != numValues; ++p) columns[p + 3] = values[p];
        outfile.writeRow(columns);
    }
    outfile.close();
}

////////////////////////////////////////////////////////////////////
