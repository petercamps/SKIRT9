/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InstrumentSystem.hpp"
#include "EscapeInstrument.hpp"
#include "FatalError.hpp"
#include "Log.hpp"

////////////////////////////////////////////////////////////////////

void InstrumentSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();

    // separate the escape instruments and arrange the peel-off instruments in groups with the same sight line
    size_t numPeelOff = 0;
    for (Instrument* instrument : _instruments)
    {
        if (auto escape = dynamic_cast<EscapeInstrument*>(instrument))
        {
            _escapeInstruments.push_back(escape);
        }
        else if (auto peelOff = dynamic_cast<PeelOffInstrument*>(instrument))
        {
            numPeelOff++;
            auto group = std::find_if(_peelOffGroups.begin(), _peelOffGroups.end(),
                                      [peelOff](const vector<PeelOffInstrument*>& candidate) {
                                          return peelOff->hasSameSightLine(candidate.front());
                                      });
            if (group != _peelOffGroups.end())
                group->push_back(peelOff);
            else
                _peelOffGroups.push_back({peelOff});
        }
        else
        {
            throw FATALERROR("Instrument is neither a peel-off instrument nor an escape instrument");
        }
    }

    // log the number of groups if this differs from the number of peel-off instruments
    if (_peelOffGroups.size() < numPeelOff)
        find<Log>()->info("The " + std::to_string(numPeelOff) + " instruments share peel-off photon packets along "
                          + std::to_string(_peelOffGroups.size()) + " sight lines");
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
