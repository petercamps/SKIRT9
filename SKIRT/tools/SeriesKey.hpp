/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SERIESKEY_HPP
#define SERIESKEY_HPP

#include "Basics.hpp"
#include <functional>
class SimulationItem;

////////////////////////////////////////////////////////////////////

/** A SeriesKey identifies a series in the IterationHistory. It consists of the simulation item on
    whose behalf the series is kept, usually the client itself, and an integer id chosen by the
    client, unique for that item, so that a client can keep several series. Because each item in
    the simulation hierarchy is a distinct object, two clients of the same type automatically have
    distinct keys. */
struct SeriesKey
{
    const SimulationItem* item{nullptr};  // the item on whose behalf the series is kept
    int id{0};                            // an id chosen by the client, unique for that item

    /** This function returns true if this key precedes the other key in the order of the item
        pointers and, for the same item, of the ids, so that keys can index a map. */
    bool operator<(const SeriesKey& other) const
    {
        if (item != other.item) return std::less<const SimulationItem*>()(item, other.item);
        return id < other.id;
    }
};

////////////////////////////////////////////////////////////////////

#endif
