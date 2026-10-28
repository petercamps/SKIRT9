/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SecondarySource.hpp"
#include <atomic>

////////////////////////////////////////////////////////////////////

SecondarySource::SecondarySource(SimulationItem* parent)
{
    parent->addChild(this);

    // because virtual functions don't work properly from with a constructor,
    // this calls setupSelfBefore/After in base classes but NOT in subclasses
    setup();
}

////////////////////////////////////////////////////////////////////

void SecondarySource::beginPreparation()
{
    // the counter is shared by all secondary sources in the process, and starts at one, so that zero never matches
    static std::atomic<size_t> counter{0};
    _preparation = ++counter;
}

////////////////////////////////////////////////////////////////////
