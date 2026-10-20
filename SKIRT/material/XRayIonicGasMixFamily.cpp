
/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "XRayIonicGasMixFamily.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "XRayIonicGasMix.hpp"

////////////////////////////////////////////////////////////////////

XRayIonicGasMixFamily::~XRayIonicGasMixFamily()
{
    for (XRayIonicGasMix* mix : _mixes) delete mix;
}

////////////////////////////////////////////////////////////////////

void XRayIonicGasMixFamily::setupSelfBefore()
{
    MaterialMixFamily::setupSelfBefore();

    // parse the required ion names from the ions property
    string ionString = StringUtils::squeeze(ions());
    if (ionString.empty()) throw FATALERROR("No ions specified");
    _ionNames = StringUtils::split(ionString, ",");
    for (string& ionName : _ionNames) ionName = StringUtils::squeeze(ionName);

    // convert our enum to XRayIonicGasMix's enum
    switch (electronScattering())
    {
        case ElectronScattering::None: _boundElectrons = XRayIonicGasMix::ElectronScattering::None; break;
        case ElectronScattering::Free: _boundElectrons = XRayIonicGasMix::ElectronScattering::Free; break;
        case ElectronScattering::FreeWithPolarization:
            _boundElectrons = XRayIonicGasMix::ElectronScattering::FreeWithPolarization;
            break;
        case ElectronScattering::FreeBound: _boundElectrons = XRayIonicGasMix::ElectronScattering::FreeBound; break;
    }

    // create a default empty mix, with zero temperature and abundances
    Array parameters(0., _ionNames.size() + 1);
    mix(parameters);
}

////////////////////////////////////////////////////////////////////

XRayIonicGasMixFamily::Capabilities XRayIonicGasMixFamily::capabilities() const
{
    Capabilities capabilities;
    capabilities.materialType = MaterialType::Gas;
    capabilities.polarizedScattering =
        electronScattering() == ElectronScattering::FreeWithPolarization || resonantScattering();
    capabilities.resonantScattering = resonantScattering();
    capabilities.scatteringDispersion = true;
    capabilities.scatteringEmulatesSecondaryEmission = true;
    return capabilities;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> XRayIonicGasMixFamily::parameterInfo() const
{
    vector<SnapshotParameter> descriptors;
    descriptors.push_back(SnapshotParameter::temperature());
    for (const string& ionName : _ionNames) descriptors.push_back(SnapshotParameter::custom(ionName));
    return descriptors;
}

////////////////////////////////////////////////////////////////////

const MaterialMix* XRayIonicGasMixFamily::mix(const Array& parameters)
{
    // split the parameters in the temperature and the abundances
    double T = parameters[0];
    vector<double> abundances(begin(parameters) + 1, end(parameters));

    // look for duplicates
    for (const XRayIonicGasMix* mix : _mixes)
    {
        if (mix->abundances() == abundances && mix->temperature() == T) return mix;
    }

    XRayIonicGasMix* mix = new XRayIonicGasMix(this, ions(), abundances, T, _boundElectrons, _resonantScattering, true);
    _mixes.push_back(mix);

    return mix;
}

////////////////////////////////////////////////////////////////////

const MaterialMix* XRayIonicGasMixFamily::mix()
{
    // the default mix is created during setup; the capabilities() function does not need it
    if (_mixes.empty()) throw FATALERROR("The default mix of an XRayIonicGasMixFamily is requested before setup");
    return _mixes[0];
}

////////////////////////////////////////////////////////////////////
