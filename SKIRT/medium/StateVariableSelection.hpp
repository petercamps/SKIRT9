/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STATEVARIABLESELECTION_HPP
#define STATEVARIABLESELECTION_HPP

#include "IterationHistory.hpp"
class MediumSystem;

////////////////////////////////////////////////////////////////////

/** A StateVariableSelection instance selects a custom medium state variable by its short name on
    behalf of a dynamic grid refinement criterion or condition (see the RefinementCriterion and
    RefinementCondition classes), and provides the value of that variable in each spatial cell.

    The variable is taken from the first medium component whose material mix offers a scalar custom
    state variable with the given name. The selection also determines the iteration loops in which
    the value of the variable changes: the primary and merged loops for a component with a primary
    dynamic medium state, the secondary and merged loops for a component with a secondary dynamic
    medium state, and, if the simulation has dynamic medium state recipes, the primary and merged
    loops for any component, because the recipes may update the state of any component in these
    loops. */
class StateVariableSelection
{
public:
    /** This function selects the variable with the specified name for the specified simulation
        item, which must be part of the medium system. It must be called during setup of the
        spatial grid or later, when the media have been set up. If no medium component offers a
        scalar custom state variable with the specified name, the function throws a fatal error
        listing the names that are available. If \em mustChange is true, the function also throws
        a fatal error if the value of the variable does not change in any iteration loop. */
    void select(const SimulationItem* item, const string& name, bool mustChange);

    /** This function returns the value of the selected variable in the spatial cell with index
        \f$m\f$. */
    double value(int m) const;

    /** This function returns the number density of the medium component that provides the
        selected variable in the spatial cell with index \f$m\f$. */
    double numberDensity(int m) const;

    /** This function returns true if the value of the selected variable changes in the specified
        iteration loop, and false otherwise. */
    bool changesIn(IterationHistory::Loop loop) const;

    /** This function returns a description of the selected variable, including the medium
        component that provides it. */
    string description() const;

private:
    const MediumSystem* _ms{nullptr};
    string _name;
    int _h{-1};
    int _index{-1};
    bool _primary{false};
    bool _secondary{false};
};

////////////////////////////////////////////////////////////////////

#endif
