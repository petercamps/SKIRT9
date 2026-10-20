/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDMEDIUMVELOCITYPROBE_HPP
#define IMPORTEDMEDIUMVELOCITYPROBE_HPP

#include "InputModelFormProbe.hpp"

////////////////////////////////////////////////////////////////////

/** ImportedMediumVelocityProbe probes the bulk velocity of each imported medium component in the
    simulation for which that information is available. The probe uses the data as represented by
    the imported snapshot, without involving the spatial grid of the simulation. The velocity is
    (mass or number) density-weighted where necessary, i.e. when the probe is associated with a
    form that projects the quantity along a path, or when two or more smoothed particles in the
    imported data overlap. */
class ImportedMediumVelocityProbe : public InputModelFormProbe
{
    ITEM_CONCRETE(ImportedMediumVelocityProbe, InputModelFormProbe, "imported medium: bulk velocity")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ImportedMediumVelocityProbe, "ImportedMedium&MediumVelocity")
    ITEM_END()

    //======================== Other Functions =======================

protected:
    /** This function probes the specified imported medium component. */
    void probeImportedMedium(const string& sh, const ImportedMedium* medium, const Snapshot* snapshot) override;
};

////////////////////////////////////////////////////////////////////

#endif
