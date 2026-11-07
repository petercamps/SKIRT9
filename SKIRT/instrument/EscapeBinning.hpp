/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ESCAPEBINNING_HPP
#define ESCAPEBINNING_HPP

#include "Direction.hpp"
#include "SimulationItem.hpp"
class ColumnOutFile;

////////////////////////////////////////////////////////////////////

/** EscapeBinning is an abstract class representing a binning of the escape directions for an
    EscapeInstrument. A binning divides the sky, i.e. the set of all directions in the model
    coordinate system, into a number of bins, each covering a part of the sky. The inclination
    \f$\theta\f$ of a direction is its angle with the z axis, and its azimuth \f$\varphi\f$ is
    measured from the x axis in the xy plane, so that a photon packet escaping in a direction that
    falls in a given bin reaches a distant instrument whose inclination and azimuth fall in that
    bin.

    Subclasses implement the functions that determine the bin of a direction and that describe
    each bin: its sky fraction, a human-readable description for the metadata of output files,
    and the values that define it, for the bin table written by the escape instrument. A subclass
    may also arrange its bins in a two-dimensional frame for the data cubes written by the escape
    instrument, by overriding both the frameLayout() and framePixel() functions; by default, the
    bins are arranged in a single row. */
class EscapeBinning : public SimulationItem
{
    ITEM_ABSTRACT(EscapeBinning, SimulationItem, "a binning of the escape directions")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the number of bins. */
    virtual int numBins() const = 0;

    /** This function returns the index of the bin containing the specified direction, or -1 if
        the direction is not recorded. The function must be thread-safe, i.e. it may not change
        the state of the receiving object. */
    virtual int bin(Direction bfk) const = 0;

    /** This function returns the fraction of the sky covered by the bin with the specified index,
        i.e. the solid angle of the bin divided by \f$4\pi\f$. */
    virtual double skyFraction(int bin) const = 0;

    /** This function returns a human-readable description of the bin with the specified index,
        for use in the metadata of output files, for example "inclination 45.6 deg to 53.1 deg". */
    virtual string describeBin(int bin) const = 0;

    /** This function adds the columns that list the values defining each bin, such as its range
        of angles, to the specified bin table, which already has columns for the bin index and the
        sky fraction. The default implementation adds no columns. */
    virtual void addBinTableColumns(ColumnOutFile& file) const;

    /** This function returns the values defining the bin with the specified index, in the order
        of the columns added by the addBinTableColumns() function. The default implementation
        returns an empty list. */
    virtual vector<double> binTableValues(int bin) const;

    /** This structure describes the frame in which the bins are arranged for the data cubes
        written by the escape instrument: the number of pixels in each direction, the increments
        between subsequent pixels and the center of the frame along each axis, and the name of the
        SKIRT quantity for these axis values, or the empty string for dimensionless values. */
    struct FrameLayout
    {
        int numPixelsX{0};
        int numPixelsY{0};
        double incrementX{1.};
        double incrementY{1.};
        double centerX{0.};
        double centerY{0.};
        string quantityXY;
    };

    /** This function returns the layout of the frame for the data cubes. The default
        implementation describes a single row of pixels, one for each bin, with the bin index as
        dimensionless axis value. */
    virtual FrameLayout frameLayout() const;

    /** This function returns the index of the pixel in the frame for the data cubes corresponding
        to the bin with the specified index, with the index along the x-axis varying most rapidly.
        The default implementation returns the bin index. */
    virtual int framePixel(int bin) const;
};

////////////////////////////////////////////////////////////////////

#endif
