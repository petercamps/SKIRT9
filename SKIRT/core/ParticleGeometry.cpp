/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParticleGeometry.hpp"
#include "ParticleSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* ParticleGeometry::createAndOpenSnapshot()
{
    // create and open the snapshot
    auto snapshot = new ParticleSnapshot;
    snapshot->open(this, filename(), "smoothed particles");

    // configure the mass column (position and size columns are configured by the snapshot itself)
    snapshot->importMass();

    // set the smoothing kernel
    snapshot->setSmoothingKernel(smoothingKernel());
    return snapshot;
}

////////////////////////////////////////////////////////////////////
