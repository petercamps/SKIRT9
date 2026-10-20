/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedMediumMetallicityProbe.hpp"
#include "EntityCollection.hpp"
#include "ImportedMedium.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedMediumMetallicityProbe::probeImportedMedium(const string& sh, const ImportedMedium* /*medium*/,
                                                         const Snapshot* snapshot)
{
    if (snapshot->hasMetallicity())
    {
        // define call-back functions to retrieve the probed quantity and corresponding weight for a given entity
        auto getValue = [](const Snapshot* snapshot, int m) { return snapshot->metallicity(m); };
        auto getWeight = [](const Snapshot* snapshot, int m) { return snapshot->density(m); };

        // construct a bridge and produce output
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity(sh + "_Z", "dimensionless", "metallicity", "density-weighted metallicity",
                             vector<const Snapshot*>{snapshot}, getValue, getWeight);
    }
}

////////////////////////////////////////////////////////////////////
