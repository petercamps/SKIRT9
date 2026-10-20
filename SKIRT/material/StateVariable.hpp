/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STATEVARIABLE_HPP
#define STATEVARIABLE_HPP

#include "Basics.hpp"
#include "QuantityKind.hpp"

////////////////////////////////////////////////////////////////////

/** StateVariable is a helper class for holding metadata about a medium state variable. For
    example, the MaterialMix::specificStateVariableInfo() function returns a list of StateVariable
    instances (see the StateVariableList class) describing the specific state variables used by the
    receiving material mix. This allows the MediumSystem class to allocate storage for the
    appropriate set of state variables, and it allows probing the relevant medium state variables
    for output.

    For each spatial cell in the simulation, the medium state includes a set of \em common state
    variables shared by all medium components and a set of \em specific state variables for each
    individual medium component. Refer to the MediumState class for a list of state variables that
    may be present in each of these sets.

    Each StateVariable instance includes an identifier specifying one of the supported state
    variables. Multiple custom variables are distinguished by an index \f$0 \le k < K\f$, where
    \f$K\f$ is the total number of custom values; the StateVariableList class assigns these indices
    when the variables are added to the list of a material mix. The StateVariable instance also
    includes a human-readable description of the quantity being represented, a unit
    system quantity definition (string), and a text output format specifier. This information
    allows custom variables to be output by generic probes or to be handled by portions of the code
    that cooperate specifically with this type of material mix, such as recipes for adjusting the
    medium state based on the calculated radiation field.

    A custom variable can also have a short name and a kind. The name consists of letters, digits,
    and underscores, and must be unique among the variables of a material mix. Two names are
    standard, so that the rest of the code can find these quantities in any material mix that
    stores them: \c T for the temperature of the material, and \c Z for its metallicity. The
    temperature() and metallicity() factory functions return custom variables with these names.
    The kind of a custom variable is intensive unless it has been marked as categorical or
    extensive (see the QuantityKind enumeration).

    A custom variable can also represent an array of values rather than a single value. Such a
    variable occupies the specified number of consecutive custom indices, starting at its own
    index. All elements of the array share the description, quantity, format, name, and kind of
    the variable.

    This class includes explicit support for the known state variables. Use one of the factory
    functions to construct a StateVariable instance. */
class StateVariable
{
    // ================== Data types ==================

public:
    /** This enumeration lists the identifiers for the supported state variables as indicated in
        the table in the class header. */
    enum class Identifier { Volume, BulkVelocity, MagneticField, NumberDensity, Custom };

    // ================== Constructing ==================

private:
    /** This constructor initializes a state variiable object with the given information for each
        field. The constructor is private; to construct a StateVariable instance, use one of
        the factory functions instead. */
    StateVariable(Identifier identifier, int customIndex, string description, string quantity, char format);

public:
    /** This function returns a StateVariable instance of type Volume. */
    static StateVariable volume();

    /** This function returns a StateVariable instance of type BulkVelocity. */
    static StateVariable bulkVelocity();

    /** This function returns a StateVariable instance of type MagneticField. */
    static StateVariable magneticField();

    /** This function returns a StateVariable instance of type NumberDensity. */
    static StateVariable numberDensity();

    /** This function returns a StateVariable instance of type Custom that holds the metallicity of
        the material, with the standard name \c Z, the description "metallicity", and a
        dimensionless quantity. */
    static StateVariable metallicity();

    /** This function returns a StateVariable instance of type Custom that holds the temperature of
        the material, with the standard name \c T, the description "temperature", and the quantity
        "temperature". */
    static StateVariable temperature();

    /** This function returns a StateVariable instance of type Custom with the specified descriptive
        information. The latter includes a human-readable description of the quantity being
        represented, a unit system quantity definition, and a text output format specifier, i.e.
        'd' for integer notation (even if the value is stored as a double), 'f' for fixed point
        notation, 'e' for scientific notation, and 'g' for the most concise 'f' or 'e'. The default
        format specifier is 'e'. The variable receives its custom index when it is added to a
        StateVariableList. */
    static StateVariable custom(const string& description, const string& quantity, char format = 'e');

    /** This function returns a copy of this custom variable with the specified short name. The
        name must consist of one or more letters, digits, and underscores. The function throws a
        fatal error if this is not the case, or if this variable is not of type Custom. */
    StateVariable named(const string& name) const;

    /** This function returns a copy of this custom variable marked as categorical. The function
        throws a fatal error if this variable is not of type Custom. */
    StateVariable categorical() const;

    /** This function returns a copy of this custom variable marked as extensive. The function
        throws a fatal error if this variable is not of type Custom. */
    StateVariable extensive() const;

    /** This function returns a copy of this custom variable that represents an array with the
        specified number of elements, which may be zero. The function throws a fatal error if this
        variable is not of type Custom, or if the number of elements is negative. */
    StateVariable array(int count) const;

    // ================== Querying ==================

public:
    /** This function returns the identifier for the state variable. All custom variables have the
        same identifier; they can be told apart through the customIndex() function(). */
    Identifier identifier() const { return _identifier; }

    /** This function returns the custom index for the state variable if it is of type Custom and
        it has been added to a StateVariableList, or zero otherwise. For an array variable, this
        is the index of its first element. */
    int customIndex() const { return _customIndex; }

    /** This function returns human-readable description for the state variable. */
    const string& description() const { return _description; }

    /** This function returns the string defining the physical quantity represented by the state
        variable. */
    const string& quantity() const { return _quantity; }

    /** This function returns the format specifier for the state variable, i.e. 'd' for integer
        notation (even if the value is stored as a double), 'f' for fixed point notation, 'e' for
        scientific notation, and 'g' for the most concise 'f' or 'e'. */
    char format() const { return _format; }

    /** This function returns the short name of the state variable, or the empty string if it has
        no name. */
    const string& name() const { return _name; }

    /** This function returns the kind of the state variable, which is QuantityKind::Intensive
        unless the variable has been marked as categorical or extensive. */
    QuantityKind kind() const { return _kind; }

    /** This function returns true if the state variable represents an array, and false if it
        represents a single value. */
    bool isArray() const { return _isArray; }

    /** This function returns the number of custom values occupied by the state variable, i.e. the
        number of elements for an array variable, and one for any other variable. */
    int count() const { return _count; }

    /** This function returns a description of the element with the specified index of an array
        variable, i.e. the description of the variable followed by a space and the index, or the
        description of the variable if it does not represent an array. */
    string elementDescription(int k) const;

    /** This function returns true if all fields of this state variable descriptor are equal to
        those of the other descriptor, and false otherwise. */
    bool operator==(const StateVariable& other) const
    {
        return _identifier == other._identifier && _customIndex == other._customIndex
               && _description == other._description && _quantity == other._quantity && _format == other._format
               && _name == other._name && _kind == other._kind && _isArray == other._isArray && _count == other._count;
    }

    /** This function returns true if some field of this state variable descriptor differs from the
        one of the other descriptor, and false otherwise. */
    bool operator!=(const StateVariable& other) const { return !(*this == other); }

    // ================== Data members ==================

private:
    // the StateVariableList class assigns the custom index
    friend class StateVariableList;

    Identifier _identifier;
    int _customIndex;
    string _description;
    string _quantity;
    char _format;
    string _name;
    QuantityKind _kind{QuantityKind::Intensive};
    bool _isArray{false};
    int _count{1};
};

////////////////////////////////////////////////////////////////////

#endif
