/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MATERIALMIXFAMILY_HPP
#define MATERIALMIXFAMILY_HPP

#include "MaterialMixProvider.hpp"

//////////////////////////////////////////////////////////////////////

/** MaterialMixFamily is an abstract class for representing a family of material mixes with the
    same capabilities, and in particular of the same material type (i.e. either dust, electrons or
    gas; never a mixture of these), and with the same specific state variables. Each subclass
    implements a material mix family from which a particular material mix can be chosen based on
    one or more parameters. This base class offers a generic interface for obtaining information on
    the number and type of parameters appropriate for the subclass, and for actually selecting a
    material mix given specific parameter values.

    Each subclass also implements the capabilities() function declared in the MaterialMixProvider
    base class, which returns the capabilities shared by all members of the family. It derives
    them from the configuration of the family, without constructing a member, because the function
    may be called before the family has been set up. When the medium system obtains the material mix
    for each spatial cell, it verifies that each of these mixes has the capabilities returned by
    this function and the same specific state variables as the default mix returned by the mix()
    function without arguments.

    The mix(const Array&) function declared in the MaterialMixProvider base class returns the
    material mix corresponding to the given values of the imported quantities, or some default
    material mix if one or more of the values are out of range. The number and type of the values
    must match the information returned by the parameterInfo() function; if not the behavior is
    undefined. The mix() function without arguments returns a default member of the family. In both
    cases, the family retains ownership of the returned material mix, and guarantees that it will
    not be destroyed until the family itself is destroyed.

    Because each material mix is represented by a distinct C++ object, a material mix family must
    by necessity consist of a discrete number of family members. As a result, each subclass must
    implement a mechanism to map specified parameter values from their continuous ranges to one of
    a discrete set of material mixes.

    Also, because a material mix object often consumes a significant amount of memory for holding
    material properties and precalculated values, the number of material mixes that can be offered
    by a material mix family is limited by the available memory resources. */
class MaterialMixFamily : public MaterialMixProvider
{
    ITEM_ABSTRACT(MaterialMixFamily, MaterialMixProvider, "a family of material mixes")
        ATTRIBUTE_TYPE_ALLOWED_IF(MaterialMixFamily, "!NonIdentitySpatialCellLibrary")
        ATTRIBUTE_TYPE_DISPLAYED_IF(MaterialMixFamily, "Level2")
    ITEM_END()

public:
    /** This function returns true because the imported quantities select a member of the family
        for each position. */
    bool hasVariableMix() const override { return true; }
};

////////////////////////////////////////////////////////////////////

#endif
