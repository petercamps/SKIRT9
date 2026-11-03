/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustDestructionRecipe.hpp"
#include "Configuration.hpp"
#include "FragmentDustMixDecorator.hpp"
#include "Log.hpp"
#include "MaterialState.hpp"
#include "MediumSystem.hpp"

////////////////////////////////////////////////////////////////////

void DustDestructionRecipe::beginUpdate(int /*numCells*/)
{
    // store radiation field wavelength grid and the medium system
    _dlambdav = find<Configuration>()->radiationFieldWLG()->dlambdav();
    _ms = find<MediumSystem>();

    // store a pointer for each medium component: either the corresponding fragmented dust mix or a null pointer
    bool hasValidFragDecorator = false;
    bool hasInvalidFragDecorator = false;
    for (auto medium : find<MediumSystem>()->media())
    {
        auto fd = dynamic_cast<const FragmentDustMixDecorator*>(medium->mix());
        if (fd && !fd->hasDynamicDensities())
        {
            fd = nullptr;
            hasInvalidFragDecorator = true;
        }
        if (fd) hasValidFragDecorator = true;
        _fdv.push_back(fd);
    }

    // inform user of problems
    if (!hasValidFragDecorator)
    {
        if (hasInvalidFragDecorator)
            find<Log>()->error("Cannot destroy dust for fragmented dust mix without dynamic densities");
        else
            find<Log>()->error("Cannot destroy dust without fragmented dust mix");
    }
    else if (hasInvalidFragDecorator)
        find<Log>()->warning("Cannot destroy dust for fragmented dust mix without dynamic densities");
}

////////////////////////////////////////////////////////////////////

UpdateStatus DustDestructionRecipe::update(MaterialState* state, const Array& Jv)
{
    UpdateStatus status;

    // act only if this component has a fragmented dust mix and if this cell contains dust
    auto fd = _fdv[state->mediumIndex()];
    if (fd && state->numberDensity() > 0.)
    {
        // get the relative noise of the radiation field in this cell, if requested
        double sigma = significance() > 0. ? _ms->radiationFieldNoise(state->cellIndex()) : 0.;

        // loop over all fragments
        int numFrags = fd->numPopulations();
        for (int f = 0; f != numFrags; ++f)
        {
            // act only if this fragment contains dust
            if (state->custom(f) > 0.)
            {
                // get the grain type and radius for this fragment
                bool graphite = fd->populationIsGraphite(f);
                double a = fd->populationGrainRadius(f);

                // get the equilibrium temperature for this fragment in the cell's radiation field
                double T = fd->populationTemperature(f, Jv);

                // call the destruction recipe to determine the density fraction
                double fraction = densityFraction(graphite, a, Jv, T);

                // if the new fraction differs by a small tolerance from the previous one, update the medium state
                double difference = abs(fraction - state->custom(numFrags + f));
                if (difference > 1e-6)  // should equal minimum value of densityFractionTolerance()
                {
                    state->setCustom(numFrags + f, fraction);

                    // if the difference is smaller than the configured tolerance, widened by the change that the
                    // noise of the radiation field could cause, consider it to be converged
                    double allowance = sigma > 0. ? noiseAllowance(fd, f, Jv, sigma) : 0.;
                    if (difference > densityFractionTolerance() + allowance)
                        status.updateNotConverged();
                    else
                        status.updateConverged();
                }
            }
        }
    }
    return status;
}

////////////////////////////////////////////////////////////////////

bool DustDestructionRecipe::usesNoise() const
{
    return significance() > 0.;
}

////////////////////////////////////////////////////////////////////

double DustDestructionRecipe::noiseAllowance(const FragmentDustMixDecorator* fd, int f, const Array& Jv,
                                             double sigma) const
{
    // a cell without contributions to one of the noise batches has infinite noise
    if (std::isinf(sigma)) return 1.;

    // calculate the fraction for the radiation field scaled down and up by the significance times the noise
    bool graphite = fd->populationIsGraphite(f);
    double a = fd->populationGrainRadius(f);
    double k = significance() * sigma;
    Array Jlo = Jv * max(0., 1. - k);
    Array Jhi = Jv * (1. + k);
    double flo = densityFraction(graphite, a, Jlo, fd->populationTemperature(f, Jlo));
    double fhi = densityFraction(graphite, a, Jhi, fd->populationTemperature(f, Jhi));
    return abs(fhi - flo);
}

////////////////////////////////////////////////////////////////////
