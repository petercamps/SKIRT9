/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MATERIALMIXFAMILY_HPP
#define MATERIALMIXFAMILY_HPP

#include "Array.hpp"
#include "MaterialMixProvider.hpp"
class MaterialMix;

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
    ITEM_END()

public:
    /** This function returns (a pointer to) the material mix corresponding to the given parameter
        values, or some default material mix if one or more of the parameter values are out of
        range. The material mix family retains ownership of the returned material mix, and
        guarantees that it will not be destroyed until the family itself is destroyed.

        The number and type of the custom parameters must match the information returned by the
        parameterInfo() function; if not the behavior is undefined. If metallicity and/or
        temperature are not being imported, this is indicated by a negative value. */
    virtual const MaterialMix* mix(double metallicity, double temperature, const Array& parameters) = 0;

    /** This function returns (a pointer to) a default material mix for this family. The material
        mix family retains ownership of the returned material mix, and guarantees that it will not
        be destroyed until the family itself is destroyed. */
    virtual const MaterialMix* mix() = 0;
};

////////////////////////////////////////////////////////////////////

#endif
