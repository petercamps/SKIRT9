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

StateVariable StateVariable::metallicity(int customIndex)
{
    return custom(customIndex, "metallicity", "").named("Z");
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::temperature(int customIndex)
{
    return custom(customIndex, "temperature", "temperature").named("T");
}

//////////////////////////////////////////////////////////////////////

StateVariable StateVariable::custom(int customIndex, const string& description, const string& quantity, char format)
{
    return StateVariable(Identifier::Custom, customIndex, description, quantity, format);
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
