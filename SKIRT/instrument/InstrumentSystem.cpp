/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InstrumentSystem.hpp"
#include "Log.hpp"

////////////////////////////////////////////////////////////////////

void InstrumentSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // arrange the instruments in groups with the same sight line
    for (Instrument* instrument : _instruments)
    {
        auto group = std::find_if(_sightLineGroups.begin(), _sightLineGroups.end(),
                                  [instrument](const vector<Instrument*>& candidate) {
                                      return instrument->hasSameSightLine(candidate.front());
                                  });
        if (group != _sightLineGroups.end())
            group->push_back(instrument);
        else
            _sightLineGroups.push_back({instrument});
    }

    // log the number of groups if this differs from the number of instruments
    if (_sightLineGroups.size() < _instruments.size())
        find<Log>()->info("The " + std::to_string(_instruments.size())
                          + " instruments share peel-off photon packets along "
                          + std::to_string(_sightLineGroups.size()) + " sight lines");
}

////////////////////////////////////////////////////////////////////

void InstrumentSystem::flush()
{
    for (Instrument* instrument : _instruments) instrument->flush();
}

////////////////////////////////////////////////////////////////////

void InstrumentSystem::write()
{
    for (Instrument* instrument : _instruments) instrument->write();
}

////////////////////////////////////////////////////////////////////
