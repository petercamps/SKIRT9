/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedMediumVelocityProbe.hpp"
#include "EntityCollection.hpp"
#include "ImportedMedium.hpp"
#include "ProbeFormBridge.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedMediumVelocityProbe::probeImportedMedium(const string& sh, const ImportedMedium* /*medium*/,
                                                      const Snapshot* snapshot)
{
    if (snapshot->hasVelocity())
    {
        // define call-back functions to retrieve the probed quantity and corresponding weight for a given entity
        auto getValue = [](const Snapshot* snapshot, int m) { return snapshot->velocity(m); };
        auto getWeight = [](const Snapshot* snapshot, int m) { return snapshot->density(m); };

        // construct a bridge and produce output
        ProbeFormBridge bridge(this, form());
        bridge.writeQuantity(sh + "_v", "velocity", "velocity", "density-weighted velocity",
                             vector<const Snapshot*>{snapshot}, getValue, getWeight);
    }
}

////////////////////////////////////////////////////////////////////
