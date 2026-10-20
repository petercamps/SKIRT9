/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STATEVARIABLELIST_HPP
#define STATEVARIABLELIST_HPP

#include "StateArray.hpp"
#include "StateScalar.hpp"
#include "StateVariable.hpp"

////////////////////////////////////////////////////////////////////

/** A StateVariableList instance assembles the specific state variables of a material mix. The list
    starts with the number density. Each custom variable that is added receives the next available
    custom index, or the next range of consecutive indices for an array variable, and the function
    adding the variable returns a handle for it (see the StateScalar and StateArray classes).

    A material mix keeps the list and the handles as data members. It fills them during setup, and
    returns the list from its MaterialMix::specificStateVariableInfo() function. Other code finds a
    variable through the list of a material mix by its name, as in \c index("T"), and loops over
    the variables of the list, as in

        for (const StateVariable& variable : mix->specificStateVariableInfo()) { ... }
    */
class StateVariableList
{
public:
    /** The constructor creates a list that holds just the number density. */
    StateVariableList();

    /** This function adds the specified scalar custom variable to the list, assigning it the next
        available custom index, and returns a handle for it. The function throws a fatal error if
        the variable is not a scalar custom variable, or if its name is already in use in the list.
        */
    StateScalar add(const StateVariable& variable);

    /** This function adds the specified array-valued custom variable to the list, assigning it the
        next available range of custom indices, and returns a handle for it. The function throws a
        fatal error if the variable is not an array-valued custom variable, or if its name is
        already in use in the list. */
    StateArray addArray(const StateVariable& variable);

    /** This function returns the custom index of the variable with the specified name, or -1 if
        the list has no such variable. For an array variable, this is the index of its first
        element. */
    int index(const string& name) const;

    /** This function returns the number of variables in the list, including the number density.
        An array variable counts as a single variable. */
    size_t size() const { return _variables.size(); }

    /** This function returns the total number of custom values occupied by the variables in the
        list, i.e. the number of custom indices assigned so far. */
    int numCustomValues() const { return _numCustomValues; }

    /** This function returns an iterator to the first variable in the list. */
    vector<StateVariable>::const_iterator begin() const { return _variables.cbegin(); }

    /** This function returns an iterator just beyond the last variable in the list. */
    vector<StateVariable>::const_iterator end() const { return _variables.cend(); }

    /** This function returns true if the two lists hold the same variables in the same order, and
        false otherwise. */
    bool operator==(const StateVariableList& other) const { return _variables == other._variables; }

    /** This function returns true if the two lists differ, and false otherwise. */
    bool operator!=(const StateVariableList& other) const { return !(*this == other); }

private:
    // adds the specified custom variable with the next available custom index, after verifying its name
    int addCustom(const StateVariable& variable);

    vector<StateVariable> _variables;
    int _numCustomValues{0};
};

////////////////////////////////////////////////////////////////////

#endif
