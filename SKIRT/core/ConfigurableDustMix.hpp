/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CONFIGURABLEDUSTMIX_HPP
#define CONFIGURABLEDUSTMIX_HPP

#include "MultiGrainDustMix.hpp"
#include "GrainPopulation.hpp"

////////////////////////////////////////////////////////////////////

/** The ConfigurableDustMix class represents a fully user-configurable dust mix described by one or
    more dust grain populations. Specifically, the class can be configured with a list of
    GrainPopulation instances, each of which represents a particular dust grain population with
    configurable grain composition, grain size distribution, and size bin discretization. During
    setup, this list of grain populations is simply passed to the base class. */
class ConfigurableDustMix : public MultiGrainDustMix
{
    ITEM_CONCRETE(ConfigurableDustMix, MultiGrainDustMix, "a configurable dust mix with one or more grain populations")

    PROPERTY_ITEM_LIST(populations, GrainPopulation, "the grain populations")
        ATTRIBUTE_DEFAULT_VALUE(populations, "GrainPopulation")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function adds the configured grain populations to the dust mix. */
    void setupSelfBefore() override;
};

////////////////////////////////////////////////////////////////////

#endif