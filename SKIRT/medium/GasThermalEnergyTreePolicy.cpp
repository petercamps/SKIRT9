/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GasThermalEnergyTreePolicy.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "Medium.hpp"
#include "MediumSystem.hpp"
#include "TreeNodeEvaluation.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // the number of positions used to estimate the average temperature of a medium component
    const int numTemperatureSamples = 100000;
}

////////////////////////////////////////////////////////////////////

void GasThermalEnergyTreePolicy::setupSelfBefore()
{
    TreePolicy::setupSelfBefore();

    // select the gas medium components that import a temperature; don't setup the medium system because we are
    // part of it
    string ignored;
    if (auto ms = find<MediumSystem>(false))
    {
        for (int h : ms->mediumIndices(MaterialMix::MaterialType::Gas))
        {
            Medium* medium = ms->media()[h];
            if (medium->hasTemperature())
            {
                _hv.push_back(h);
                _media.push_back(medium);
            }
            else
                ignored += (ignored.empty() ? "" : ", ") + std::to_string(h);
        }
    }
    if (_hv.empty()) throw FATALERROR(type() + " requires a gas medium component that imports a temperature");
    if (!ignored.empty())
        find<Log>()->warning(type()
                             + " ignores the gas medium components that do not import a temperature: " + ignored);

    // estimate the average temperature and the total thermal energy of these components (without the constant factor)
    double total = 0.;
    for (Medium* medium : _media)
    {
        double sumT = 0.;
        for (int i = 0; i != numTemperatureSamples; ++i) sumT += medium->temperature(medium->generatePosition());
        _averageTv.push_back(sumT / numTemperatureSamples);
        total += medium->number() * _averageTv.back();
    }
    if (!(total > 0.)) throw FATALERROR("The total thermal energy of the gas is not positive");
    _maxEnergy = maxFraction() * total;
}

////////////////////////////////////////////////////////////////////

bool GasThermalEnergyTreePolicy::needsSubdivide(TreeNodeEvaluation& node) const
{
    const auto& positionv = node.samplePositions();
    int numSamples = positionv.size();
    double energy = 0.;
    for (size_t k = 0; k != _hv.size(); ++k)
    {
        // get the number of entities in the node, which is calculated exactly if possible
        double number = node.density(_hv[k]) * node.box().volume();
        if (number > 0.)
        {
            // determine the average temperature in the node, weighted by the sampled density, or use the average
            // temperature of the medium component if none of the samples hits any material
            const auto& densityv = node.componentDensitySamples(_hv[k]);
            double sumn = 0.;
            double sumnT = 0.;
            for (int i = 0; i != numSamples; ++i)
            {
                if (densityv[i] > 0.)
                {
                    sumn += densityv[i];
                    sumnT += densityv[i] * _media[k]->temperature(positionv[i]);
                }
            }
            energy += number * (sumn > 0. ? sumnT / sumn : _averageTv[k]);
        }
    }
    return energy > _maxEnergy;
}

////////////////////////////////////////////////////////////////////
