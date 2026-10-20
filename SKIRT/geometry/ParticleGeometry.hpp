/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARTICLEGEOMETRY_HPP
#define PARTICLEGEOMETRY_HPP

#include "ImportedGeometry.hpp"
#include "SmoothingKernel.hpp"

////////////////////////////////////////////////////////////////////

/** A ParticleGeometry instance represents a 3D geometry with a spatial density distribution
    described by a list of smoothed particles. The particle data is usually extracted from a
    cosmological simulation snapshot, and it must be provided in a column text file formatted as
    described below. The total mass in the geometry is normalized to unity after importing the
    data.

    Refer to the description of the TextInFile class for information on overall formatting and on
    how to include header lines specifying the units for each column in the input file. In case the
    input file has no unit specifications, the default units mentioned below are used instead. The
    input file should contain 5 columns:

    \f[ x\,(\mathrm{pc}) \quad y\,(\mathrm{pc}) \quad z\,(\mathrm{pc}) \quad h\,(\mathrm{pc}) \quad
    M\,(\mathrm{M}_\odot) \f]

    The first three columns are the \f$x\f$, \f$y\f$ and \f$z\f$ coordinates of the particle, the
    fourth column is the particle smoothing length \f$h\f$, and the fifth column is the mass
    \f$M\f$ of the particle. The mass units are actually irrelevant because the total mass in the
    geometry will be normalized to unity after importing the data. However, the import procedure
    still insists on knowing the units. */
class ParticleGeometry : public ImportedGeometry
{
    ITEM_CONCRETE(ParticleGeometry, ImportedGeometry, "a geometry imported from smoothed particle data")

        PROPERTY_ITEM(smoothingKernel, SmoothingKernel, "the kernel for interpolating the smoothed particles")
        ATTRIBUTE_DEFAULT_VALUE(smoothingKernel, "CubicSplineSmoothingKernel")
        ATTRIBUTE_DISPLAYED_IF(smoothingKernel, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new ParticleSnapshot object, calls its open() function,
        configures it to import a mass column, passes the smoothing kernel selected by the user to
        it, and finally returns a pointer to the object. Ownership of the Snapshot object is
        transferred to the caller. */
    Snapshot* createAndOpenSnapshot() override;
};

////////////////////////////////////////////////////////////////////

#endif
