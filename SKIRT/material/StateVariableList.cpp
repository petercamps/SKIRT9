/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StateVariableList.hpp"
#include "FatalError.hpp"

//////////////////////////////////////////////////////////////////////

StateVariableList::StateVariableList() : _variables{StateVariable::numberDensity()} {}

//////////////////////////////////////////////////////////////////////

StateScalar StateVariableList::add(const StateVariable& variable)
{
    if (variable.identifier() != StateVariable::Identifier::Custom || variable.isArray())
        throw FATALERROR("Only a scalar custom state variable can be added with StateVariableList::add()");
    return StateScalar(addCustom(variable));
}

//////////////////////////////////////////////////////////////////////

StateArray StateVariableList::addArray(const StateVariable& variable)
{
    if (variable.identifier() != StateVariable::Identifier::Custom || !variable.isArray())
        throw FATALERROR("Only an array-valued custom state variable can be added with StateVariableList::addArray()");
    return StateArray(addCustom(variable), variable.count());
}

//////////////////////////////////////////////////////////////////////

int StateVariableList::index(const string& name) const
{
    for (const StateVariable& variable : _variables)
        if (!name.empty() && variable.name() == name) return variable.customIndex();
    return -1;
}

//////////////////////////////////////////////////////////////////////

int StateVariableList::addCustom(const StateVariable& variable)
{
    if (index(variable.name()) >= 0)
        throw FATALERROR("The state variable name '" + variable.name() + "' is used more than once");

    int customIndex = _numCustomValues;
    _variables.push_back(variable);
    _variables.back()._customIndex = customIndex;
    _numCustomValues += variable.count();
    return customIndex;
}

//////////////////////////////////////////////////////////////////////
