/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OpticalDepthTreePolicy.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"
#include "TreeNodeEvaluation.hpp"

////////////////////////////////////////////////////////////////////

void OpticalDepthTreePolicy::setupSelfBefore()
{
    TreePolicy::setupSelfBefore();

    // determine the indices of the medium components to be included
    // (don't setup the medium system because we are part of it)
    auto ms = find<MediumSystem>(false);
    int numMedia = ms ? ms->media().size() : 0;
    if (materialType() == MaterialType::All)
    {
        for (int h = 0; h != numMedia; ++h) _hv.push_back(h);
        if (_hv.empty()) throw FATALERROR(type() + " for all media requires at least one medium component");
    }
    else
    {
        auto mixType = MaterialMix::MaterialType::Dust;
        if (materialType() == MaterialType::Electrons) mixType = MaterialMix::MaterialType::Electrons;
        if (materialType() == MaterialType::Gas) mixType = MaterialMix::MaterialType::Gas;
        if (ms) _hv = ms->mediumIndices(mixType);
        if (_hv.empty())
            throw FATALERROR(type() + " for " + TreeNodeEvaluation::materialTypeName(mixType)
                             + " requires a medium component of that material type");
    }

    // determine the extinction coefficient of each medium component: per unit mass for dust, per entity otherwise
    for (int h : _hv)
    {
        auto mix = ms->media()[h]->mix();
        double sigma = mix->sectionExt(wavelength());
        _kappav.push_back(mix->materialType() == MaterialMix::MaterialType::Dust ? sigma / mix->mass() : sigma);
    }
}

////////////////////////////////////////////////////////////////////

bool OpticalDepthTreePolicy::needsSubdivide(TreeNodeEvaluation& node) const
{
    // add the contributions of the medium components to the optical depth
    double tau = 0.;
    int numComponents = _hv.size();
    for (int k = 0; k != numComponents; ++k) tau += _kappav[k] * node.density(_hv[k]) * node.box().diagonal();
    return tau > maxOpticalDepth();
}

////////////////////////////////////////////////////////////////////

Range OpticalDepthTreePolicy::wavelengthRange() const
{
    return Range(wavelength(), wavelength());
}

////////////////////////////////////////////////////////////////////
