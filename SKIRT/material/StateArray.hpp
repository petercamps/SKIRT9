/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STATEARRAY_HPP
#define STATEARRAY_HPP

#include "MaterialState.hpp"

////////////////////////////////////////////////////////////////////

/** A StateArray instance is a handle that refers to an array-valued custom state variable of a
    material mix, i.e. a variable that occupies a number of consecutive custom indices. The
    StateVariableList class hands out the handle when the variable is added to the list of the
    material mix. The material mix keeps the handle as a data member, and applies it to a material
    state to obtain a view on the elements of the variable in the corresponding spatial cell and
    medium component:

        auto kappa = _opacityExt(state);
        double k = kappa[j] + (kappa[j + 1] - kappa[j]) * fraction;

    The view refers directly to the values stored in the medium state, so that no values are
    copied, and the position of the variable is calculated only once for all elements. A view
    obtained from a constant material state is read-only. A default constructed handle refers to
    no variable, and must not be applied to a material state. */
class StateArray
{
public:
    /** The default constructor creates a handle that refers to no variable. */
    StateArray() {}

    /** This function returns a read-only view on the elements of the variable in the specified
        material state. */
    Span<const double> operator()(const MaterialState* state) const { return state->customs(_index, _count); }

    /** This function returns a writable view on the elements of the variable in the specified
        material state. */
    Span<double> operator()(MaterialState* state) const { return state->customs(_index, _count); }

    /** This function returns the custom index of the first element of the variable, or -1 if the
        handle refers to no variable. */
    int index() const { return _index; }

    /** This function returns the number of elements of the variable. */
    int count() const { return _count; }

private:
    // only the StateVariableList class creates handles that refer to a variable
    friend class StateVariableList;

    /** This constructor creates a handle that refers to the variable with the specified custom
        index and number of elements. */
    StateArray(int index, int count) : _index{index}, _count{count} {}

    int _index{-1};
    int _count{0};
};

////////////////////////////////////////////////////////////////////

#endif
