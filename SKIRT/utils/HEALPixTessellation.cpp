/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "HEALPixTessellation.hpp"

////////////////////////////////////////////////////////////////////

HEALPixTessellation::HEALPixTessellation(int order) : _order{order}, _Nside{1 << order} {}

////////////////////////////////////////////////////////////////////

int HEALPixTessellation::numPixelsInRing(int j) const
{
    int r = j + 1;  // one-based ring number
    if (r < _Nside) return 4 * r;
    if (r > 3 * _Nside) return 4 * (4 * _Nside - r);
    return 4 * _Nside;
}

////////////////////////////////////////////////////////////////////

std::pair<int, int> HEALPixTessellation::locate(double theta, double phi) const
{
    // the code below was mostly copied from healpix_base.cc in the official HEALPix repository
    // it corresponds to the loc2pix function using the RING scheme
    // unlike the loc2pix function, we do not compute a full RING pixel index, but restrict
    // ourselves to a computation of the ring index and pixel-in-ring index, which we respectively
    // use as vertical and horizontal pixel index
    double z = cos(theta);
    double za = abs(z);
    double tt = fmod(2. * phi / M_PI, 4.);

    // i: pixel-in-ring index (horizontal pixel index)
    // j: ring index (vertical pixel index)
    int i, j;
    // first figure out if we are in the equatorial or polar region
    if (za <= 2. / 3.)
    {
        // equatorial region: all rings have 4*_Nside pixels
        double t1 = _Nside * (0.5 + tt);
        double t2 = 0.75 * _Nside * z;
        int jp = static_cast<int>(floor(t1 - t2));
        int jm = static_cast<int>(floor(t1 + t2));

        // we first compute the ring index as the offset in the equatorial region (j in [1, 2*_Nside+1])
        j = _Nside + 1 + jp - jm;
        int kshift = 1 - (j & 1);
        int temp = jp + jm + kshift + 1 + 7 * _Nside;
        i = (_order > 0) ? ((temp >> 1) & (4 * _Nside - 1)) : ((temp >> 1) % (4 * _Nside));
        // we now add the correct offset to the ring index to account for the rings in the north polar region
        j += _Nside - 2;
    }
    else
    {
        // polar region: the number of pixels per ring depends on how far the ring is from the pole
        double tp = tt - static_cast<int>(tt);
        double tmp = (za < 0.99) ? (_Nside * sqrt(3. * (1. - za))) : (_Nside * sin(theta) / sqrt((1. + za) / 3.));

        int jp = static_cast<int>(tp * tmp);
        int jm = static_cast<int>((1. - tp) * tmp);

        // we first compute the ring index as an offset from the pole (j in [1, _Nside-1])
        j = jp + jm + 1;
        i = static_cast<int>(tt * j);
        // we now convert the ring index to the actual ring index, and distinguish between north and south pole
        if (z < 0)
        {
            j = 4 * _Nside - j - 1;
        }
        else
        {
            --j;
        }
    }
    return {i, j};
}

////////////////////////////////////////////////////////////////////

int HEALPixTessellation::ringPixelIndex(int i, int j) const
{
    int r = j + 1;  // one-based ring number
    if (r < _Nside) return 2 * r * (r - 1) + i;
    if (r <= 3 * _Nside) return 2 * _Nside * (_Nside - 1) + (r - _Nside) * 4 * _Nside + i;
    int rs = 4 * _Nside - r;  // ring number counted from the south pole
    return numPixels() - 2 * rs * (rs + 1) + i;
}

////////////////////////////////////////////////////////////////////

std::pair<int, int> HEALPixTessellation::ringPixelLocation(int p) const
{
    int ncap = 2 * _Nside * (_Nside - 1);  // number of pixels in the north polar cap
    int r = 0;                             // one-based ring number
    int i = 0;                             // pixel-in-ring index
    if (p < ncap)
    {
        // north polar cap: ring r starts at pixel 2r(r-1)
        r = static_cast<int>(0.5 * (1. + sqrt(1. + 2. * p)));
        while (2 * r * (r - 1) > p) --r;
        while (2 * (r + 1) * r <= p) ++r;
        i = p - 2 * r * (r - 1);
    }
    else if (p < numPixels() - ncap)
    {
        // equatorial region: all rings have 4*_Nside pixels
        int q = p - ncap;
        r = _Nside + q / (4 * _Nside);
        i = q % (4 * _Nside);
    }
    else
    {
        // south polar cap: ring rs counted from the south pole starts at pixel numPixels - 2rs(rs+1)
        int ps = numPixels() - p;  // in [1, ncap]
        int rs = static_cast<int>(0.5 * (-1. + sqrt(1. + 2. * ps)));
        while (2 * rs * (rs + 1) < ps) ++rs;
        while (rs > 1 && 2 * (rs - 1) * rs >= ps) --rs;
        r = 4 * _Nside - rs;
        i = p - (numPixels() - 2 * rs * (rs + 1));
    }
    return {i, r - 1};
}

////////////////////////////////////////////////////////////////////

std::pair<double, double> HEALPixTessellation::pixelCenter(int i, int j) const
{
    int r = j + 1;  // one-based ring number
    double z = 0.;
    double phi = 0.;
    if (r < _Nside)
    {
        // north polar cap
        z = 1. - static_cast<double>(r) * r / (3. * _Nside * _Nside);
        phi = (i + 0.5) * M_PI / (2. * r);
    }
    else if (r > 3 * _Nside)
    {
        // south polar cap
        int rs = 4 * _Nside - r;
        z = -1. + static_cast<double>(rs) * rs / (3. * _Nside * _Nside);
        phi = (i + 0.5) * M_PI / (2. * rs);
    }
    else
    {
        // equatorial region: alternate rings are shifted by half a pixel
        z = (2 * _Nside - r) * 2. / (3. * _Nside);
        double shift = ((r + _Nside) & 1) ? 0. : 0.5;
        phi = (i + shift) * M_PI / (2. * _Nside);
    }
    return {acos(z), phi};
}

////////////////////////////////////////////////////////////////////
