/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ImportedMedium.hpp"
#include "Configuration.hpp"
#include "Snapshot.hpp"

////////////////////////////////////////////////////////////////////

void ImportedMedium::setupSelfAfter()
{
    Medium::setupSelfAfter();

    // create the snapshot with preconfigured mass or density column
    _snapshot = createAndOpenSnapshot();

    // add optional standard columns if applicable
    if (hasVelocity()) _snapshot->importVelocity();
    if (_importMagneticField) _snapshot->importMagneticField();

    // add the columns for the quantities requested by the material mix or material mix family
    vector<SnapshotParameter> parameterInfo = _materialMix->parameterInfo();
    if (!parameterInfo.empty()) _snapshot->importParameters(parameterInfo);

    // set the density policy
    _snapshot->setMassDensityPolicy(_massFraction);

    // notify about building search data structures if needed
    if (find<Configuration>()->snapshotsNeedGetEntities()) _snapshot->setNeedGetEntities();

    // read the data from file
    _snapshot->readAndClose();
}

//////////////////////////////////////////////////////////////////////

ImportedMedium::~ImportedMedium()
{
    delete _snapshot;
}

//////////////////////////////////////////////////////////////////////

int ImportedMedium::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

const MaterialMix* ImportedMedium::mix(Position bfr) const
{
    if (_materialMix->hasVariableMix())
    {
        Array params;
        _snapshot->parameters(bfr, params);
        return _materialMix->mix(params);
    }
    else
        return _materialMix->mix();
}

//////////////////////////////////////////////////////////////////////

const MaterialMix* ImportedMedium::mix() const
{
    return _materialMix->mix();
}

//////////////////////////////////////////////////////////////////////

MaterialMixProvider::Capabilities ImportedMedium::capabilities() const
{
    return _materialMix->capabilities();
}

//////////////////////////////////////////////////////////////////////

bool ImportedMedium::hasVariableMix() const
{
    return _materialMix->hasVariableMix();
}

////////////////////////////////////////////////////////////////////

bool ImportedMedium::hasUniformMass() const
{
    return _materialMix->hasUniformMass();
}

//////////////////////////////////////////////////////////////////////

bool ImportedMedium::hasVelocity() const
{
    if (_importVelocity)
    {
        // refuse velocity for oligochromatic simulations
        // (this function is called from Configure so we cannot precompute this during setup)
        auto config = find<Configuration>();
        if (!config->oligochromatic()) return true;
    }
    return false;
}

//////////////////////////////////////////////////////////////////////

Vec ImportedMedium::bulkVelocity(Position bfr) const
{
    return hasVelocity() ? _snapshot->velocity(bfr) : Vec();
}

//////////////////////////////////////////////////////////////////////

bool ImportedMedium::hasMagneticField() const
{
    return _importMagneticField;
}

//////////////////////////////////////////////////////////////////////

Vec ImportedMedium::magneticField(Position bfr) const
{
    return _importMagneticField ? _snapshot->magneticField(bfr) : Vec();
}

////////////////////////////////////////////////////////////////////

bool ImportedMedium::hasParameters() const
{
    return _snapshot->hasParameters();
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> ImportedMedium::parameterInfo() const
{
    return _materialMix->parameterInfo();
}

////////////////////////////////////////////////////////////////////

void ImportedMedium::parameters(Position bfr, Array& params) const
{
    if (hasParameters())
        _snapshot->parameters(bfr, params);
    else
        params.resize(0);
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::numberDensity(Position bfr) const
{
    // with a uniform mass per entity, the default mix suffices, so that no family member is looked up
    double result = _snapshot->density(bfr);
    if (!_snapshot->holdsNumber()) result /= (hasUniformMass() ? mix() : mix(bfr))->mass();
    return result;
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::number() const
{
    double result = _snapshot->mass();
    if (!_snapshot->holdsNumber()) result /= mix()->mass();
    return result;
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::massDensity(Position bfr) const
{
    // with a uniform mass per entity, the default mix suffices, so that no family member is looked up
    double result = _snapshot->density(bfr);
    if (_snapshot->holdsNumber()) result *= (hasUniformMass() ? mix() : mix(bfr))->mass();
    return result;
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::mass() const
{
    double result = _snapshot->mass();
    if (_snapshot->holdsNumber()) result *= mix()->mass();
    return result;
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::opticalDepthX(double lambda) const
{
    double result = _snapshot->SigmaX() * mix()->sectionExt(lambda);
    if (!_snapshot->holdsNumber()) result /= mix()->mass();
    return result;
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::opticalDepthY(double lambda) const
{
    double result = _snapshot->SigmaY() * mix()->sectionExt(lambda);
    if (!_snapshot->holdsNumber()) result /= mix()->mass();
    return result;
}

////////////////////////////////////////////////////////////////////

double ImportedMedium::opticalDepthZ(double lambda) const
{
    double result = _snapshot->SigmaZ() * mix()->sectionExt(lambda);
    if (!_snapshot->holdsNumber()) result /= mix()->mass();
    return result;
}

//////////////////////////////////////////////////////////////////////

Position ImportedMedium::generatePosition() const
{
    return _snapshot->generatePosition();
}

//////////////////////////////////////////////////////////////////////

int ImportedMedium::numSites() const
{
    return _snapshot->numEntities();
}

//////////////////////////////////////////////////////////////////////

Position ImportedMedium::sitePosition(int index) const
{
    return _snapshot->position(index);
}

//////////////////////////////////////////////////////////////////////

const Snapshot* ImportedMedium::snapshot() const
{
    return _snapshot;
}

//////////////////////////////////////////////////////////////////////
