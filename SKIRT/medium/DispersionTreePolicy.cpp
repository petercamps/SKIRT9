/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DispersionTreePolicy.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"
#include "TreeNodeEvaluation.hpp"

////////////////////////////////////////////////////////////////////

void DispersionTreePolicy::setupSelfBefore()
{
    TreePolicy::setupSelfBefore();

    // determine the selected material type
    switch (materialType())
    {
        case MaterialType::Dust: _type = MaterialMix::MaterialType::Dust; break;
        case MaterialType::Electrons: _type = MaterialMix::MaterialType::Electrons; break;
        case MaterialType::Gas: _type = MaterialMix::MaterialType::Gas; break;
    }

    // verify that there are medium components of that type (don't setup the medium system because we are part of it)
    auto ms = find<MediumSystem>(false);
    if (!ms || ms->mediumIndices(_type).empty())
        throw FATALERROR(type() + " for " + TreeNodeEvaluation::materialTypeName(_type)
                         + " requires a medium component of that material type");
}

////////////////////////////////////////////////////////////////////

bool DispersionTreePolicy::needsSubdivide(TreeNodeEvaluation& node) const
{
    // determine the smallest and largest sampled density in the node
    const auto& samples = node.densitySamples(_type);
    auto [minIt, maxIt] = std::minmax_element(samples.cbegin(), samples.cend());
    double rhomin = *minIt;
    double rhomax = *maxIt;

    double q = rhomax > 0 ? (rhomax - rhomin) / rhomax : 0.;
    return q > maxDispersion();
}

////////////////////////////////////////////////////////////////////
