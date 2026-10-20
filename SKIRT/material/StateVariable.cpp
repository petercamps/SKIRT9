/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StateVariable.hpp"
#include "FatalError.hpp"
#include <cctype>

//////////////////////////////////////////////////////////////////////

StateVariable::StateVariable(Identifier identifier, int customIndex, string description, string quantity, char format)
    : _identifier(identifier), _customIndex(customIndex), _description(description), _quantity(quantity),
      _format(format)
{}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::volume()
{
    return StateVariable(Identifier::Volume, 0, "volume", "volume", 'e');
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::bulkVelocity()
{
    return StateVariable(Identifier::BulkVelocity, 0, "bulk velocity", "velocity", 'e');
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::magneticField()
{
    return StateVariable(Identifier::MagneticField, 0, "magnetic field", "magneticfield", 'e');
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::numberDensity()
{
    return StateVariable(Identifier::NumberDensity, 0, "number density", "numbervolumedensity", 'e');
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::metallicity()
{
    return custom("metallicity", "").named("Z");
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::temperature()
{
    return custom("temperature", "temperature").named("T");
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::custom(const string& description, const string& quantity, char format)
{
    return StateVariable(Identifier::Custom, 0, description, quantity, format);
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::named(const string& name) const
{
    if (_identifier != Identifier::Custom) throw FATALERROR("Only a custom state variable can have a name");
    bool valid = !name.empty();
    for (char c : name)
        if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_') valid = false;
    if (!valid) throw FATALERROR("Invalid state variable name '" + name + "'");

    StateVariable result = *this;
    result._name = name;
    return result;
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::categorical() const
{
    if (_identifier != Identifier::Custom) throw FATALERROR("Only a custom state variable can be categorical");
    StateVariable result = *this;
    result._kind = QuantityKind::Categorical;
    return result;
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::extensive() const
{
    if (_identifier != Identifier::Custom) throw FATALERROR("Only a custom state variable can be extensive");
    StateVariable result = *this;
    result._kind = QuantityKind::Extensive;
    return result;
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::array(int count) const
{
    if (_identifier != Identifier::Custom) throw FATALERROR("Only a custom state variable can be an array");
    if (count < 0) throw FATALERROR("The number of elements of an array state variable cannot be negative");
    StateVariable result = *this;
    result._isArray = true;
    result._count = count;
    return result;
}

//////////////////////////////////////////////////////////////////////

string StateVariable::elementDescription(int k) const
{
    return _isArray ? _description + " " + std::to_string(k) : _description;
}

//////////////////////////////////////////////////////////////////////
