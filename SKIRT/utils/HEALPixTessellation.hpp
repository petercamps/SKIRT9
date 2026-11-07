/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HEALPIXTESSELLATION_HPP
#define HEALPIXTESSELLATION_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** HEALPixTessellation is a helper class that locates directions on a HEALPix tessellation of the
    sphere, which divides the sphere into pixels of equal area. It is used by the
    HEALPixSkyInstrument and HEALPixEscapeBinning classes.

    The class is based on the official HEALPix algorithm by Górski, K. M. et al. (2005)
    (https://ui.adsabs.harvard.edu/abs/2005ApJ...622..759G/abstract), as available from
    https://healpix.sourceforge.io, but only uses the relevant part of the algorithm.

    The tessellation of order \f$k\f$ has resolution parameter \f$N_\mathrm{side}=2^k\f$ and
    contains \f$12N_\mathrm{side}^2\f$ pixels, with centers on \f$4N_\mathrm{side}-1\f$ rings of
    equal latitude. A ring with (zero-based) ring index \f$j\f$ corresponds to the ring number
    \f$j+1\f$ in the HEALPix papers, counted from the north pole (\f$\theta=0\f$). The rings with
    \f$N_\mathrm{side}-1\le j\le 3N_\mathrm{side}-1\f$ contain \f$4N_\mathrm{side}\f$ pixels each;
    the number of pixels in the polar rings decreases by 4 for every ring closer to the pole.
    Within a ring, the (zero-based) pixel-in-ring index \f$i\f$ increases with the azimuth
    \f$\varphi\f$.

    The class offers the two indices \f$(i,j)\f$ of the pixel containing a given direction, which
    can be used as the horizontal and vertical pixel coordinates in an image of
    \f$4N_\mathrm{side}\times(4N_\mathrm{side}-1)\f$ pixels, with unused pixels in the polar
    regions. It also offers the conversion of these indices to the pixel index in the HEALPix \em
    ring ordering scheme, and the direction of the pixel center. */
class HEALPixTessellation
{
public:
    /** The constructor creates a tessellation of the specified order, which must be in the range
        from 0 to 15. */
    explicit HEALPixTessellation(int order = 0);

    /** This function returns the order of the tessellation. */
    int order() const { return _order; }

    /** This function returns the resolution parameter \f$N_\mathrm{side}=2^k\f$, i.e. the number
        of pixels along each side of a HEALPix base pixel. */
    int nside() const { return _Nside; }

    /** This function returns the number of pixels in the tessellation,
        \f$12N_\mathrm{side}^2\f$. */
    int numPixels() const { return 12 * _Nside * _Nside; }

    /** This function returns the number of rings in the tessellation, \f$4N_\mathrm{side}-1\f$,
        which is also the number of rows in the image layout. */
    int numRings() const { return 4 * _Nside - 1; }

    /** This function returns the largest number of pixels in a ring, \f$4N_\mathrm{side}\f$, which
        is also the number of columns in the image layout. */
    int maxPixelsInRing() const { return 4 * _Nside; }

    /** This function returns the number of pixels in the ring with the specified (zero-based)
        ring index. */
    int numPixelsInRing(int j) const;

    /** This function returns the solid angle of each pixel, \f$\pi/(3N_\mathrm{side}^2)\f$. */
    double pixelSolidAngle() const { return M_PI / (3. * _Nside * _Nside); }

    /** This function determines the pixel containing the direction with the specified polar angle
        \f$\theta\in[0,\pi]\f$ and azimuth angle \f$\varphi\in[0,2\pi]\f$, and returns its
        pixel-in-ring index \f$i\f$ and ring index \f$j\f$, in that order. */
    std::pair<int, int> locate(double theta, double phi) const;

    /** This function returns the index of the pixel with the specified pixel-in-ring index and ring
        index in the HEALPix ring ordering scheme. */
    int ringPixelIndex(int i, int j) const;

    /** This function returns the pixel-in-ring index and the ring index, in that order, of the
        pixel with the specified index in the HEALPix ring ordering scheme. */
    std::pair<int, int> ringPixelLocation(int p) const;

    /** This function returns the polar angle \f$\theta\f$ and the azimuth angle \f$\varphi\f$,
        in that order, of the center of the pixel with the specified pixel-in-ring index and ring
        index. */
    std::pair<double, double> pixelCenter(int i, int j) const;

private:
    int _order{0};
    int _Nside{1};
};

////////////////////////////////////////////////////////////////////

#endif
