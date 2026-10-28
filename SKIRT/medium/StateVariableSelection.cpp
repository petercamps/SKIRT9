/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StateVariableSelection.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "MediumSystem.hpp"
#include "StateVariableList.hpp"
#include "StringUtils.hpp"

////////////////////////////////////////////////////////////////////

void StateVariableSelection::select(const SimulationItem* item, const string& name, bool mustChange)
{
    // the medium system is not set up, because the item is part of it; its media have been set up, however
    _ms = item->find<MediumSystem>(false);
    _name = name;

    // find the first medium component that offers a scalar custom variable with the specified name,
    // collecting the available names in case there is none
    vector<string> available;
    int numMedia = _ms->media().size();
    for (int h = 0; h != numMedia && _h < 0; ++h)
    {
        for (const StateVariable& variable : _ms->media()[h]->mix()->specificStateVariableInfo())
        {
            if (variable.name().empty()) continue;
            if (variable.name() == name && variable.count() == 1)
            {
                _h = h;
                _index = variable.customIndex();
                break;
            }
            available.push_back(variable.name());
        }
    }
    if (_h < 0)
    {
        string message = "No medium component offers a scalar state variable with name '" + name + "'";
        if (available.empty()) throw FATALERROR(message + "; there are no named state variables");
        throw FATALERROR(message + "; available names: " + StringUtils::join(available, ", "));
    }

    // determine the iteration loops in which the value of the variable changes
    auto config = item->find<Configuration>();
    bool merged = config->hasMergedIterations();
    switch (_ms->media()[_h]->capabilities().dynamicMediumState)
    {
        case MaterialMix::DynamicStateType::None: break;
        case MaterialMix::DynamicStateType::Primary: _primary = true; break;
        case MaterialMix::DynamicStateType::Secondary: _secondary = true; break;
        case MaterialMix::DynamicStateType::PrimaryIfMergedIterations:
            if (merged)
                _primary = true;
            else
                _secondary = true;
            break;
    }
    if (config->hasDynamicStateRecipes()) _primary = true;

    if (mustChange && !_primary && !_secondary)
        throw FATALERROR("The state variable '" + name + "' of medium component " + std::to_string(_h)
                         + " does not change between iterations, because the component has no dynamic medium state");
}

////////////////////////////////////////////////////////////////////

double StateVariableSelection::value(int m) const
{
    return _ms->custom(m, _h, _index);
}

////////////////////////////////////////////////////////////////////

double StateVariableSelection::numberDensity(int m) const
{
    return _ms->numberDensity(m, _h);
}

////////////////////////////////////////////////////////////////////

bool StateVariableSelection::changesIn(IterationHistory::Loop loop) const
{
    switch (loop)
    {
        case IterationHistory::Loop::None: return false;
        case IterationHistory::Loop::Primary: return _primary;
        case IterationHistory::Loop::Secondary: return _secondary;
        case IterationHistory::Loop::Merged: return _primary || _secondary;
    }
    return false;
}

////////////////////////////////////////////////////////////////////

string StateVariableSelection::description() const
{
    return _name + " of medium component " + std::to_string(_h);
}

////////////////////////////////////////////////////////////////////
