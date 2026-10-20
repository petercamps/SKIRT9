/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MaterialMix.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

void MaterialMix::setupSelfBefore()
{
    MaterialMixProvider::setupSelfBefore();

    _random = find<Random>();
    _config = find<Configuration>();
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasPolarizedScattering() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasPolarizedAbsorption() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasPolarizedEmission() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasResonantScattering() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasNegativeExtinction() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasStochasticDustEmission() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasExtraSpecificState() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasScatteringDispersion() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::scatteringEmulatesSecondaryEmission() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

MaterialMix::DynamicStateType MaterialMix::hasDynamicMediumState() const
{
    return DynamicStateType::None;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasContinuumEmission() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasLineEmission() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

MaterialMix::Capabilities MaterialMix::capabilities() const
{
    Capabilities capabilities;
    capabilities.materialType = materialType();
    capabilities.polarizedScattering = hasPolarizedScattering();
    capabilities.polarizedAbsorption = hasPolarizedAbsorption();
    capabilities.polarizedEmission = hasPolarizedEmission();
    capabilities.resonantScattering = hasResonantScattering();
    capabilities.negativeExtinction = hasNegativeExtinction();
    capabilities.stochasticDustEmission = hasStochasticDustEmission();
    capabilities.extraSpecificState = hasExtraSpecificState();
    capabilities.scatteringDispersion = hasScatteringDispersion();
    capabilities.scatteringEmulatesSecondaryEmission = scatteringEmulatesSecondaryEmission();
    capabilities.continuumEmission = hasContinuumEmission();
    capabilities.lineEmission = hasLineEmission();
    capabilities.dynamicMediumState = hasDynamicMediumState();
    return capabilities;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> MaterialMix::parameterInfo() const
{
    return vector<SnapshotParameter>();
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::hasVariableMix() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

const MaterialMix* MaterialMix::mix(const Array& /*parameters*/)
{
    return this;
}

////////////////////////////////////////////////////////////////////

const MaterialMix* MaterialMix::mix()
{
    return this;
}

////////////////////////////////////////////////////////////////////

void MaterialMix::initializeSpecificState(MaterialState* /*state*/, const Array& /*params*/) const {}

////////////////////////////////////////////////////////////////////

UpdateStatus MaterialMix::updateSpecificState(MaterialState* /*state*/, const Array& /*Jv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::isSpecificStateConverged(int /*numCells*/, int /*numUpdated*/, int /*numNotConverged*/,
                                           MaterialState* /*currentAggregate*/,
                                           MaterialState* /*previousAggregate*/) const
{
    return true;
}

////////////////////////////////////////////////////////////////////

double MaterialMix::asymmpar(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

bool MaterialMix::peeloffScattering(double& /*I*/, double& /*Q*/, double& /*U*/, double& /*V*/, double& /*lambda*/,
                                    Direction /*bfkobs*/, Direction /*bfky*/, const MaterialState* /*state*/,
                                    const PhotonPacket* /*pp*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

void MaterialMix::performScattering(double /*lambda*/, const MaterialState* /*state*/, PhotonPacket* /*pp*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

DisjointWavelengthGrid* MaterialMix::emissionWavelengthGrid() const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

Array MaterialMix::emissivity(const Array& /*Jv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

Array MaterialMix::emissionSpectrum(const MaterialState* /*state*/, const Array& /*Jv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

const Array& MaterialMix::thetaGrid() const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

const Array& MaterialMix::sectionsAbs(double /*lambda*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

const Array& MaterialMix::sectionsAbspol(double /*lambda*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

Array MaterialMix::lineEmissionCenters() const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

Array MaterialMix::lineEmissionMasses() const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

Array MaterialMix::lineEmissionSpectrum(const MaterialState* /*state*/, const Array& /*Jv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////

double MaterialMix::indicativeTemperature(const MaterialState* /*state*/, const Array& /*Jv*/) const
{
    throw FATALERROR("This function implementation should never be called");
}

////////////////////////////////////////////////////////////////////
