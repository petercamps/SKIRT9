/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MATERIALMIXPROVIDER_HPP
#define MATERIALMIXPROVIDER_HPP

#include "Array.hpp"
#include "SimulationItem.hpp"
#include "SnapshotParameter.hpp"
class MaterialMix;

////////////////////////////////////////////////////////////////////

/** MaterialMixProvider is the abstract base class of the MaterialMix and MaterialMixFamily
    classes. A material mix defines the material properties of a medium component, and a material
    mix family provides a material mix for each position, selected by quantities imported from a
    snapshot.

    This base class defines the enumerations for the material type and for the type of dynamic
    medium state, and the Capabilities structure, which describes the capabilities of a material
    mix or, for a family, of all its members. Because a family describes its capabilities without
    constructing a member, the configuration of a simulation can be established from the
    capabilities of each medium component before the medium components have been set up.

    The enumerations remain accessible through the MaterialMix class, so that code referring to,
    for example, MaterialMix::MaterialType::Dust is not affected by their move to this class. */
class MaterialMixProvider : public SimulationItem
{
    ITEM_ABSTRACT(MaterialMixProvider, SimulationItem, "a material mix or a family of material mixes")
    ITEM_END()

    //======== Data types =======

public:
    /** This enumeration lists the fundamental material types supported by the MaterialMix class
        hierarchy. */
    enum class MaterialType { Dust, Electrons, Gas };

    /** This enumeration is used to indicate whether a dynamic medium state (DMS) update algorithm
        is provided, and if so, when it should be executed:

        - \c None: no DMS update algorithm is provided.

        - \c Primary: the DMS update affects the opacity for primary emission (PDMS).

        - \c Secondary: the DMS update affects the secondary emission only (SDMS).

        - \c PrimaryIfMergedIterations: \c Primary if the simulation has merged iterations (\em
        iterateSecondaryEmission and \em includePrimaryEmission are both true), and \c Secondary
        otherwise. */
    enum class DynamicStateType { None, Primary, Secondary, PrimaryIfMergedIterations };

    /** This structure describes the capabilities of a material mix, or of all members of a
        material mix family. Each data member corresponds to the MaterialMix function with a
        similar name, as listed below; refer to the documentation of that function for its
        meaning. The default values are those returned by the default implementations of these
        functions in the MaterialMix class; the default material type, which that class does not
        provide, is arbitrary.

        Data member | MaterialMix function
        ------------|---------------------
        materialType | materialType()
        polarizedScattering | hasPolarizedScattering()
        polarizedAbsorption | hasPolarizedAbsorption()
        polarizedEmission | hasPolarizedEmission()
        resonantScattering | hasResonantScattering()
        negativeExtinction | hasNegativeExtinction()
        stochasticDustEmission | hasStochasticDustEmission()
        extraSpecificState | hasExtraSpecificState()
        scatteringDispersion | hasScatteringDispersion()
        scatteringEmulatesSecondaryEmission | scatteringEmulatesSecondaryEmission()
        continuumEmission | hasContinuumEmission()
        lineEmission | hasLineEmission()
        dynamicMediumState | hasDynamicMediumState()
        */
    struct Capabilities
    {
        MaterialType materialType{MaterialType::Dust};
        bool polarizedScattering{false};
        bool polarizedAbsorption{false};
        bool polarizedEmission{false};
        bool resonantScattering{false};
        bool negativeExtinction{false};
        bool stochasticDustEmission{false};
        bool extraSpecificState{false};
        bool scatteringDispersion{false};
        bool scatteringEmulatesSecondaryEmission{false};
        bool continuumEmission{false};
        bool lineEmission{false};
        DynamicStateType dynamicMediumState{DynamicStateType::None};

        /** This function returns true if all data members of this structure are equal to those of
            the other structure, and false otherwise. */
        bool operator==(const Capabilities& other) const;

        /** This function returns true if some data member of this structure differs from the one
            of the other structure, and false otherwise. */
        bool operator!=(const Capabilities& other) const { return !(*this == other); }
    };

    //======== Interface =======

public:
    /** This function returns the number and type of the quantities to be imported from a
        snapshot for this material mix or material mix family, as a list of SnapshotParameter
        objects. Each of these objects specifies unit information and a human-readable description
        for the quantity. For a material mix, the quantities are used to initialize the medium
        state of each spatial cell; for a family, they select a member of the family at each
        position. */
    virtual vector<SnapshotParameter> parameterInfo() const = 0;

    /** This function returns the capabilities of this material mix or, for a material mix family,
        the capabilities shared by all its members. The function may be called before the receiving
        item has been set up, and it does not construct any material mix. */
    virtual Capabilities capabilities() const = 0;

    /** This function returns true if the imported quantities select a material mix for each
        position, as for a material mix family, and false if they initialize the medium state of a
        single material mix. */
    virtual bool hasVariableMix() const = 0;

    /** This function returns (a pointer to) the material mix for the specified values of the
        imported quantities, in the order returned by the parameterInfo() function. A material mix
        returns itself, and a material mix family returns one of its members. The receiving item
        retains ownership of the returned material mix. */
    virtual const MaterialMix* mix(const Array& parameters) = 0;

    /** This function returns (a pointer to) a material mix representative of this material mix or
        family. A material mix returns itself, and a material mix family returns a default member.
        The receiving item retains ownership of the returned material mix. */
    virtual const MaterialMix* mix() = 0;
};

////////////////////////////////////////////////////////////////////

#endif
