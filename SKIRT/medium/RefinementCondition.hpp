/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef REFINEMENTCONDITION_HPP
#define REFINEMENTCONDITION_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** RefinementCondition is the abstract base class for the conditions that restrict a dynamic grid
    refinement criterion to part of the medium, for example to partly ionized gas (see the
    RefinementCriterion and DynamicRefinementOptions classes). A condition provides the current
    value of a quantity in each spatial cell, and tells whether a value averaged over the
    iterations of the averaging window of the refinement lies in its range. A criterion ignores
    every cell in which one of its conditions does not hold. A condition keeps no per-cell data of
    its own; the medium system accumulates its values in the iteration history. */
class RefinementCondition : public SimulationItem
{
    ITEM_ABSTRACT(RefinementCondition, SimulationItem, "a dynamic grid refinement condition")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function returns the current value of the quantity of the condition in the spatial
        cell with index \f$m\f$. The function may be called from multiple threads at the same time.
        */
    virtual double value(int m) const = 0;

    /** This function returns true if the specified value, averaged over the iterations of the
        averaging window, lies in the range of the condition, and false otherwise. */
    virtual bool holds(double value) const = 0;

    /** This function returns a description of the quantity of the condition, for use in log
        messages and in the description of the corresponding cell window in the iteration history.
        */
    virtual string description() const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
