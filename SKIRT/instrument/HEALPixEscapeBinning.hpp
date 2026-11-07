/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HEALPIXESCAPEBINNING_HPP
#define HEALPIXESCAPEBINNING_HPP

#include "EscapeBinning.hpp"
#include "HEALPixTessellation.hpp"

////////////////////////////////////////////////////////////////////

/** HEALPixEscapeBinning is an escape binning that divides the sky into the
    \f$12N_\mathrm{side}^2\f$ pixels of equal area of a HEALPix tessellation, with
    \f$N_\mathrm{side}=2^k\f$ for the specified order \f$k\f$. The tessellation is the same as that
    of the HEALPixSkyInstrument (see the HEALPixTessellation class), aligned with the model
    coordinates: the polar angle of the tessellation is the inclination \f$\theta\f$, and its
    azimuth angle is the azimuth \f$\varphi\in[0,2\pi]\f$. The bin index is the pixel index in the
    HEALPix \em ring ordering scheme.

    For the data cubes written by the escape instrument, the bins are arranged in the same
    two-dimensional layout as the images of the HEALPixSkyInstrument, i.e. a frame of
    \f$4N_\mathrm{side}\times(4N_\mathrm{side}-1)\f$ pixels, with the pixel-in-ring index along
    the x-axis (roughly the azimuth) and the ring index along the y-axis (roughly the inclination),
    so that the frames of the cube can be viewed as maps of the sky. The pixels of the frame that
    do not correspond to a HEALPix pixel remain zero. */
class HEALPixEscapeBinning : public EscapeBinning
{
    ITEM_CONCRETE(HEALPixEscapeBinning, EscapeBinning, "the pixels of a HEALPix tessellation")

        PROPERTY_INT(order, "HEALPix order")
        ATTRIBUTE_MIN_VALUE(order, "0")
        ATTRIBUTE_MAX_VALUE(order, "10")
        ATTRIBUTE_DEFAULT_VALUE(order, "3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function creates the HEALPix tessellation of the specified order. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns the number of bins, i.e. the number of HEALPix pixels. */
    int numBins() const override;

    /** This function returns the index of the HEALPix pixel containing the specified direction, in
        the ring ordering scheme. */
    int bin(Direction bfk) const override;

    /** This function returns the fraction of the sky covered by each HEALPix pixel, which is the
        same for all bins. */
    double skyFraction(int bin) const override;

    /** This function returns a description of the bin with the specified index, listing its ring
        index, its pixel-in-ring index, and the direction of its center. */
    string describeBin(int bin) const override;

    /** This function adds columns for the ring index, the pixel-in-ring index, and the inclination
        and azimuth of the center of each bin. */
    void addBinTableColumns(ColumnOutFile& file) const override;

    /** This function returns the ring index, the pixel-in-ring index, and the inclination and
        azimuth of the center (in output units) of the bin with the specified index. */
    vector<double> binTableValues(int bin) const override;

    /** This function describes the frame used for the data cubes, which has the same layout as
        the images of the HEALPixSkyInstrument. */
    FrameLayout frameLayout() const override;

    /** This function returns the index of the pixel in the frame for the data cubes corresponding
        to the bin with the specified index. */
    int framePixel(int bin) const override;

    //======================== Data Members ========================

private:
    HEALPixTessellation _tessellation;
};

////////////////////////////////////////////////////////////////////

#endif
