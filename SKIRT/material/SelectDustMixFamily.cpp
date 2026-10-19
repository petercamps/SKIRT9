/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SelectDustMixFamily.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void SelectDustMixFamily::setupSelfAfter()
{
    MaterialMixFamily::setupSelfAfter();

    auto capabilities = _dustMixes[0]->capabilities();
    int numMixes = _dustMixes.size();
    for (int index = 1; index != numMixes; ++index)
        if (_dustMixes[index]->capabilities() != capabilities)
            throw FATALERROR("The dust mix with index " + std::to_string(index)
                             + " in the SelectDustMixFamily has other capabilities than the dust mix with index 0");
}

////////////////////////////////////////////////////////////////////

SelectDustMixFamily::Capabilities SelectDustMixFamily::capabilities() const
{
    return _dustMixes[0]->capabilities();
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> SelectDustMixFamily::parameterInfo() const
{
    return {SnapshotParameter::custom("dustmix index").categorical()};
}

////////////////////////////////////////////////////////////////////

const MaterialMix* SelectDustMixFamily::mix(double /*metallicity*/, double /*temperature*/, const Array& parameters)
{
    long numMixes = _dustMixes.size();
    long index = max(0L, min(std::lround(parameters[0]), numMixes - 1));
    return _dustMixes[index];
}

////////////////////////////////////////////////////////////////////

const MaterialMix* SelectDustMixFamily::mix()
{
    return _dustMixes[0];
}

////////////////////////////////////////////////////////////////////
