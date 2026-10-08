/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef APERTUREINSTRUMENT_HPP
#define APERTUREINSTRUMENT_HPP

#include "DistantInstrument.hpp"

////////////////////////////////////////////////////////////////////

/** ApertureInstrument is an abstract subclass of DistantInstrument that offers an extra \em radius
    property for the user to configure a circular aperture centered on the origin of the model
    coordinate system and in the plane perpendicular to the instrument's line of sight. Photon
    packets arriving from a point that parallel projects outside of this aperture are ignored. If
    the radius is zero (the default value), the instrument does not have an aperture (or,
    equivalently, the aperture radius is infinite). */
class ApertureInstrument : public DistantInstrument
{
    ITEM_ABSTRACT(ApertureInstrument, DistantInstrument,
                  "a distant instrument that offers an optional circular aperture")

        PROPERTY_DOUBLE(radius, "the radius of the circular aperture, or zero for no aperture")
        ATTRIBUTE_QUANTITY(radius, "length")
        ATTRIBUTE_MIN_VALUE(radius, "[0")
        ATTRIBUTE_DEFAULT_VALUE(radius, "0")
        ATTRIBUTE_DISPLAYED_IF(radius, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function stores some information used by the locate() function. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function determines whether the instrument would detect a peel-off photon packet
        launched towards it from the specified position. It returns a pixel index of zero if the
        position projects inside the configured radius, or if that radius is zero, and a pixel
        index of -1 otherwise. The distance is left at its default value of infinity. */
    Detection locate(Position bfr) const override;

    //======================== Data Members ========================

private:
    // data members derived from the discoverable properties during setup, used in locate()
    double _radius2{0};
    double _costheta{0};
    double _sintheta{0};
    double _cosphi{0};
    double _sinphi{0};
};

////////////////////////////////////////////////////////////////////

#endif
