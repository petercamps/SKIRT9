/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "DustLuminosityCriterion.hpp"
#include "CellWindow.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"

////////////////////////////////////////////////////////////////////

void DustLuminosityCriterion::setupSelfAfter()
{
    RefinementCriterion::setupSelfAfter();

    // the medium system is not set up, because the criterion is part of it; the lists of medium indices are available
    _ms = find<MediumSystem>(false);
    if (_ms->dustMediumIndices().empty())
        throw FATALERROR("A dust luminosity refinement criterion requires a dust medium component");
    if (!find<Configuration>()->hasDustEmission())
        throw FATALERROR("A dust luminosity refinement criterion requires dust emission");
}

////////////////////////////////////////////////////////////////////

bool DustLuminosityCriterion::changesIn(IterationHistory::Loop loop) const
{
    // the primary radiation field is recalculated in the primary and merged loops only
    if (primaryOnly()) return loop == IterationHistory::Loop::Primary || loop == IterationHistory::Loop::Merged;
    return loop != IterationHistory::Loop::None;
}

////////////////////////////////////////////////////////////////////

void DustLuminosityCriterion::prepareValues() const
{
    auto total = _ms->totalDustAbsorbedLuminosity();
    _total = primaryOnly() ? total.first : total.first + total.second;
}

////////////////////////////////////////////////////////////////////

double DustLuminosityCriterion::value(int m) const
{
    if (_total <= 0.) return 0.;
    return (primaryOnly() ? _ms->primaryDustLuminosity(m) : _ms->dustLuminosity(m)) / _total;
}

////////////////////////////////////////////////////////////////////

double DustLuminosityCriterion::numberDensity(int m) const
{
    double n = 0.;
    for (int h : _ms->dustMediumIndices()) n += _ms->numberDensity(m, h);
    return n;
}

////////////////////////////////////////////////////////////////////

string DustLuminosityCriterion::fieldDescription() const
{
    return "fraction of the dust-absorbed luminosity";
}

////////////////////////////////////////////////////////////////////

bool DustLuminosityCriterion::needsNeighbors() const
{
    return false;
}

////////////////////////////////////////////////////////////////////

double DustLuminosityCriterion::excess(int m, const CellWindow& window, const vector<int>& /*neighbors*/) const
{
    double ratio = window.mean(m) / maxFraction();
    return ratio > 1. ? ratio : 0.;
}

////////////////////////////////////////////////////////////////////
