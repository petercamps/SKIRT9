/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STATESCALAR_HPP
#define STATESCALAR_HPP

#include "MaterialState.hpp"

////////////////////////////////////////////////////////////////////

/** A StateScalar instance is a handle that refers to a scalar custom state variable of a material
    mix. The StateVariableList class hands out the handle when the variable is added to the list of
    the material mix. The material mix keeps the handle as a data member, and applies it to a
    material state to read or write the value of the variable in the corresponding spatial cell
    and medium component:

        double fraction = _hNeutralFraction(state);
        _hNeutralFraction(state) = 1.;

    The accessors are inline, so that using a handle is as efficient as passing the custom index
    of the variable to the custom() and setCustom() functions of the material state. A material
    mix that receives a constant material state can read the value but cannot change it. A default
    constructed handle refers to no variable, and must not be applied to a material state. */
class StateScalar
{
public:
    /** The default constructor creates a handle that refers to no variable. */
    StateScalar() {}

    /** This function returns the value of the variable in the specified material state. */
    double operator()(const MaterialState* state) const { return state->custom(_index); }

    /** This function returns a writable reference to the value of the variable in the specified
        material state. */
    double& operator()(MaterialState* state) const { return state->customs(_index, 1)[0]; }

    /** This function returns the custom index of the variable, or -1 if the handle refers to no
        variable. */
    int index() const { return _index; }

private:
    // only the StateVariableList class creates handles that refer to a variable
    friend class StateVariableList;

    /** This constructor creates a handle that refers to the variable with the specified custom
        index. */
    explicit StateScalar(int index) : _index{index} {}

    int _index{-1};
};

////////////////////////////////////////////////////////////////////

#endif
