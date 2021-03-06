/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "RadialVectorField.hpp"

//////////////////////////////////////////////////////////////////////

int RadialVectorField::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

Vec RadialVectorField::vector(Position bfr) const
{
    double r = bfr.norm();
    double s = r < unityRadius() ? r / unityRadius() : 1.;
    Vec u = r > 0. ? bfr / r : Vec();
    return s * u;
}

//////////////////////////////////////////////////////////////////////
