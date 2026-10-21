/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef AGGREGATIONRULE_HPP
#define AGGREGATIONRULE_HPP

#include "Basics.hpp"
class Medium;

////////////////////////////////////////////////////////////////////

/** An AggregationRule specifies the value that the medium system records in an aggregate series
    of the IterationHistory: the volume integral over all spatial cells of the custom state
    variable with the specified index, for the specified medium component. The history only stores
    and compares the medium pointer, and never dereferences it. */
struct AggregationRule
{
    const Medium* medium{nullptr};  // the medium component
    int customIndex{0};             // the index of the custom state variable

    /** This function returns true if the two rules specify the same medium component and custom
        state variable, and false otherwise. */
    bool operator==(const AggregationRule& other) const
    {
        return medium == other.medium && customIndex == other.customIndex;
    }
};

////////////////////////////////////////////////////////////////////

#endif
