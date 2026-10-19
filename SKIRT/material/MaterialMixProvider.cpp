/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MaterialMixProvider.hpp"

////////////////////////////////////////////////////////////////////

bool MaterialMixProvider::Capabilities::operator==(const Capabilities& other) const
{
    return materialType == other.materialType && polarizedScattering == other.polarizedScattering
           && polarizedAbsorption == other.polarizedAbsorption && polarizedEmission == other.polarizedEmission
           && resonantScattering == other.resonantScattering && negativeExtinction == other.negativeExtinction
           && stochasticDustEmission == other.stochasticDustEmission && extraSpecificState == other.extraSpecificState
           && scatteringDispersion == other.scatteringDispersion
           && scatteringEmulatesSecondaryEmission == other.scatteringEmulatesSecondaryEmission
           && continuumEmission == other.continuumEmission && lineEmission == other.lineEmission
           && dynamicMediumState == other.dynamicMediumState;
}

////////////////////////////////////////////////////////////////////
