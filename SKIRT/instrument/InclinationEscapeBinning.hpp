/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INCLINATIONESCAPEBINNING_HPP
#define INCLINATIONESCAPEBINNING_HPP

#include "EscapeBinning.hpp"

////////////////////////////////////////////////////////////////////

/** InclinationEscapeBinning is an escape binning that divides the sky into bins of inclination
    and, optionally, azimuth. The inclination bins are equally spaced in \f$\cos\theta\f$, so that
    all bins cover the same solid angle, and they are ordered from the positive z axis
    (\f$\theta=0\f$) towards the negative z axis. If the \em foldHemispheres property is true, the
    directions with inclinations \f$\theta\f$ and \f$180^\circ-\theta\f$ fall in the same bin, so
    that the bins are equally spaced in \f$|\cos\theta|\f$ and each bin covers both hemispheres; the
    bins are then ordered from the z axis towards the xy plane. This option suits models that are
    symmetric with respect to the xy plane. The azimuth bins have equal width in \f$\varphi\f$,
    starting at \f$\varphi=0\f$. The bin index is \f$i_\theta N_\varphi + i_\varphi\f$, where
    \f$i_\theta\f$ and \f$i_\varphi\f$ are the indices of the inclination and azimuth bins, and
    \f$N_\varphi\f$ is the number of azimuth bins. */
class InclinationEscapeBinning : public EscapeBinning
{
    ITEM_CONCRETE(InclinationEscapeBinning, EscapeBinning, "bins of inclination and, optionally, azimuth")

        PROPERTY_INT(numInclinationBins, "the number of inclination bins, equally spaced in cos(θ)")
        ATTRIBUTE_MIN_VALUE(numInclinationBins, "1")
        ATTRIBUTE_MAX_VALUE(numInclinationBins, "100000")
        ATTRIBUTE_DEFAULT_VALUE(numInclinationBins, "10")

        PROPERTY_BOOL(foldHemispheres, "combine the directions θ and 180°-θ in the same bin")
        ATTRIBUTE_DEFAULT_VALUE(foldHemispheres, "false")

        PROPERTY_INT(numAzimuthBins, "the number of azimuth bins, of equal width in φ")
        ATTRIBUTE_MIN_VALUE(numAzimuthBins, "1")
        ATTRIBUTE_MAX_VALUE(numAzimuthBins, "100000")
        ATTRIBUTE_DEFAULT_VALUE(numAzimuthBins, "1")
        ATTRIBUTE_DISPLAYED_IF(numAzimuthBins, "Dimension3")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the number of bins, i.e. the product of the numbers of inclination
        and azimuth bins. */
    int numBins() const override;

    /** This function returns the index of the bin containing the specified direction. */
    int bin(Direction bfk) const override;

    /** This function returns the fraction of the sky covered by the bin with the specified index,
        which is the same for all bins. */
    double skyFraction(int bin) const override;

    /** This function returns a description of the bin with the specified index, listing its
        inclination range and, if there is more than one azimuth bin, its azimuth range. */
    string describeBin(int bin) const override;

    /** This function adds columns for the inclination range and, if there is more than one
        azimuth bin, the azimuth range of each bin. If the hemispheres are folded, the inclination
        range lists the part of the bin in the hemisphere with positive z. */
    void addBinTableColumns(ColumnOutFile& file) const override;

    /** This function returns the inclination range and, if applicable, the azimuth range of the
        bin with the specified index, in output units. */
    vector<double> binTableValues(int bin) const override;

private:
    /** This function returns the inclination and azimuth ranges of the bin with the specified
        index, in radians, as the minimum and maximum inclination followed by the minimum and
        maximum azimuth. If the hemispheres are folded, the inclination range is the part of the
        bin in the hemisphere with positive z. */
    std::tuple<double, double, double, double> binRanges(int bin) const;
};

////////////////////////////////////////////////////////////////////

#endif
