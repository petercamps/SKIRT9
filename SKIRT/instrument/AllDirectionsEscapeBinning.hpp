/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ALLDIRECTIONSESCAPEBINNING_HPP
#define ALLDIRECTIONSESCAPEBINNING_HPP

#include "EscapeBinning.hpp"

////////////////////////////////////////////////////////////////////

/** AllDirectionsEscapeBinning is an escape binning with a single bin that covers the whole sky. An
    escape instrument with this binning records the angle-averaged spectrum of the radiation
    escaping from the model, such as the omni component of a torus model. */
class AllDirectionsEscapeBinning : public EscapeBinning
{
    ITEM_CONCRETE(AllDirectionsEscapeBinning, EscapeBinning, "all directions together in a single bin")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the number of bins, which is one. */
    int numBins() const override;

    /** This function returns the index of the single bin, zero, for any direction. */
    int bin(Direction bfk) const override;

    /** This function returns the fraction of the sky covered by the single bin, which is one. */
    double skyFraction(int bin) const override;

    /** This function returns a description of the single bin. */
    string describeBin(int bin) const override;
};

////////////////////////////////////////////////////////////////////

#endif
