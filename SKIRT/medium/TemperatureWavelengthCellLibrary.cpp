/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TemperatureWavelengthCellLibrary.hpp"
#include "Configuration.hpp"
#include "Log.hpp"
#include "Medium.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProcessManager.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"
#include "WavelengthGrid.hpp"

////////////////////////////////////////////////////////////////////

int TemperatureWavelengthCellLibrary::numEntries() const
{
    return _numTemperatures * _numWavelengths;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // a dust medium component with its weight, i.e. its total mass in the simulation relative to that of all dust
    struct DustComponent
    {
        const MaterialMix* mix{nullptr};
        double weight{0.};
    };

    // returns the dust medium components, weighted by their total mass in the simulation
    vector<DustComponent> dustComponents(const MediumSystem* ms)
    {
        int numCells = ms->numCells();
        vector<DustComponent> components;
        double totalMass = 0.;
        for (int h = 0; h != ms->numMedia(); ++h)
        {
            if (ms->isDust(h))
            {
                double mass = 0.;
                for (int m = 0; m != numCells; ++m) mass += ms->massDensity(m, h) * ms->volume(m);
                components.push_back({ms->media()[h]->mix(), mass});
                totalMass += mass;
            }
        }

        // normalize the weights to a sum of one; if the total mass is zero, give all components the same weight
        for (auto& component : components)
            component.weight = totalMass > 0. ? component.weight / totalMass : 1. / components.size();
        return components;
    }
}

////////////////////////////////////////////////////////////////////

vector<int> TemperatureWavelengthCellLibrary::mapping(const Array& bv) const
{
    // get the radiation field wavelength grid and the medium system
    auto wavelengthGrid = find<Configuration>()->radiationFieldWLG();
    auto ms = find<MediumSystem>();
    int numCells = ms->numCells();

    // get the dust components, weighted by their total mass in the simulation
    vector<DustComponent> components = dustComponents(ms);

    // tabulate the mass-weighted absorption mass coefficient of the dust components times the bin width, and times
    // the wavelength, on the radiation field wavelength grid
    int numWavelengths = wavelengthGrid->numBins();
    Array kappadlambdav(numWavelengths);
    Array kappadlambdalambdav(numWavelengths);
    for (int ell = 0; ell != numWavelengths; ++ell)
    {
        double lambda = wavelengthGrid->wavelength(ell);
        double kappa = 0.;
        for (const auto& component : components)
        {
            double mu = component.mix->mass();
            if (mu > 0.) kappa += component.weight * component.mix->sectionAbs(lambda) / mu;
        }
        kappadlambdav[ell] = kappa * wavelengthGrid->effectiveWidth(ell);
        kappadlambdalambdav[ell] = kappadlambdav[ell] * lambda;
    }

    // calculate the indicative temperature and wavelength for all spatial cells in parallel, each cell on a single
    // process, and sum the results over the processes, so that they are identical on all processes
    Array Tv(numCells);
    Array lambdav(numCells);
    find<ParallelFactory>()->parallelDistributed()->call(numCells, [&](size_t firstIndex, size_t numIndices) {
        for (size_t m = firstIndex; m != firstIndex + numIndices; ++m)
        {
            // ignore cells that won't be used by the caller
            if (bv[m])
            {
                // calculate both properties for the dust components in the radiation field of the cell; the material
                // state is not needed: a dust mix does not use it, and a cell library cannot be combined with a
                // fragmented dust mix, which does
                const Array& Jv = ms->meanIntensity(m);
                double T = 0.;
                for (const auto& component : components)
                    T += component.weight * component.mix->indicativeTemperature(nullptr, Jv);
                double sumbot = (kappadlambdav * Jv).sum();
                double lambda = sumbot > 0. ? (kappadlambdalambdav * Jv).sum() / sumbot : 0.;

                // ignore cells with meaningless property values
                if (T > 0. && lambda > 0.)
                {
                    Tv[m] = T;
                    lambdav[m] = lambda;
                }
            }
        }
    });
    ProcessManager::sumToAll(Tv);
    ProcessManager::sumToAll(lambdav);

    // track the minimum and maximum values
    double Tmin = DBL_MAX;
    double Tmax = 0.0;
    double lambdamin = DBL_MAX;
    double lambdamax = 0.0;
    for (int m = 0; m != numCells; ++m)
    {
        if (Tv[m] > 0. && lambdav[m] > 0.)
        {
            Tmin = min(Tmin, Tv[m]);
            Tmax = max(Tmax, Tv[m]);
            lambdamin = min(lambdamin, lambdav[m]);
            lambdamax = max(lambdamax, lambdav[m]);
        }
    }

    // log the property ranges
    auto log = find<Log>();
    auto units = find<Units>();
    log->info("  Indicative temperatures vary"
              " from T = "
              + StringUtils::toString(units->otemperature(Tmin), 'g', 5) + " " + units->utemperature()
              + " to T = " + StringUtils::toString(units->otemperature(Tmax), 'g', 5) + " " + units->utemperature());
    log->info("  Indicative wavelengths vary"
              " from λ = "
              + StringUtils::toString(units->owavelength(lambdamin), 'g', 5) + " " + units->uwavelength()
              + " to λ = " + StringUtils::toString(units->owavelength(lambdamax), 'g', 5) + " " + units->uwavelength());

    // determine for every dust cell m the corresponding library entry n
    double dT = (Tmax - Tmin) / _numTemperatures;
    double loglambdamin = log10(lambdamin);
    double loglambdamax = log10(lambdamax);
    double dloglambda = (loglambdamax - loglambdamin) / _numWavelengths;
    vector<int> nv(numCells);
    for (int m = 0; m != numCells; ++m)
    {
        if (Tv[m] > 0. && lambdav[m] > 0.)
        {
            // if all qualifying cells share the same temperature and/or wavelength, the corresponding
            // bin width is zero; put them all in bin 0 for that dimension
            double T = Tv[m];
            int i = dT > 0. ? max(0, min(_numTemperatures - 1, static_cast<int>((T - Tmin) / dT))) : 0;

            double lambda = lambdav[m];
            double loglambda = log10(lambda);
            int j = dloglambda > 0.
                        ? max(0, min(_numWavelengths - 1, static_cast<int>((loglambda - loglambdamin) / dloglambda)))
                        : 0;

            nv[m] = i + _numTemperatures * j;
        }
        else
        {
            nv[m] = -1;
        }
    }

    return nv;
}

////////////////////////////////////////////////////////////////////
