/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ITERATIONOPTIONS_HPP
#define ITERATIONOPTIONS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The IterationOptions class simply offers a number of options for configuring iterations during
    primary and/or secondary emission. These options are relevant only when the simulation has a
    dynamic medium state and/or a dynamic secondary emission.

    The number of photon packets launched in each primary emission iteration is the number of
    primary packets multiplied by \em primaryIterationPacketsMultiplier. If \em
    primaryIterationInitialPacketsFraction is below one, the first iteration launches that fraction
    of this number, and each next iteration \em primaryIterationPacketsRamp times as many, up to the
    full number. Because the iterations with fewer packets are noisier, the primary emission loop
    does not count as converged before the number of packets has reached the full number. If \em
    primaryIterationPacketsPerCell is nonzero, each iteration instead launches that number of
    packets per spatial cell, bounded by the number for the first iteration and the full number.
    This is meant for a grid with dynamic refinement, whose number of cells grows during the
    iterations: the noise in the field of a refinement criterion depends on the number of packets
    per cell. */
class IterationOptions : public SimulationItem
{
    ITEM_CONCRETE(IterationOptions, SimulationItem, "a set of options for primary and/or secondary emission iterations")

        PROPERTY_INT(minPrimaryIterations, "the minimum number of iterations during primary emission")
        ATTRIBUTE_MIN_VALUE(minPrimaryIterations, "1")
        ATTRIBUTE_MAX_VALUE(minPrimaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minPrimaryIterations, "1")
        ATTRIBUTE_RELEVANT_IF(minPrimaryIterations, "IteratePrimary")

        PROPERTY_INT(maxPrimaryIterations, "the maximum number of iterations during primary emission")
        ATTRIBUTE_MIN_VALUE(maxPrimaryIterations, "1")
        ATTRIBUTE_MAX_VALUE(maxPrimaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxPrimaryIterations, "10")
        ATTRIBUTE_RELEVANT_IF(maxPrimaryIterations, "IteratePrimary")

        PROPERTY_INT(minSecondaryIterations, "the minimum number of iterations during secondary emission")
        ATTRIBUTE_MIN_VALUE(minSecondaryIterations, "1")
        ATTRIBUTE_MAX_VALUE(minSecondaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minSecondaryIterations, "1")
        ATTRIBUTE_RELEVANT_IF(minSecondaryIterations, "IterateSecondary")

        PROPERTY_INT(maxSecondaryIterations, "the maximum number of iterations during secondary emission")
        ATTRIBUTE_MIN_VALUE(maxSecondaryIterations, "0")
        ATTRIBUTE_MAX_VALUE(maxSecondaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxSecondaryIterations, "10")
        ATTRIBUTE_RELEVANT_IF(maxSecondaryIterations, "IterateSecondary")

        PROPERTY_BOOL(includePrimaryEmission, "include primary emission in the secondary emission iterations")
        ATTRIBUTE_DEFAULT_VALUE(includePrimaryEmission, "false")
        ATTRIBUTE_RELEVANT_IF(includePrimaryEmission, "IterateSecondary")
        ATTRIBUTE_DISPLAYED_IF(includePrimaryEmission, "Level3")

        PROPERTY_DOUBLE(primaryIterationPacketsMultiplier,
                        "the multiplier on the number of packets for each primary iteration")
        ATTRIBUTE_MIN_VALUE(primaryIterationPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(primaryIterationPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(primaryIterationPacketsMultiplier, "1")
        ATTRIBUTE_RELEVANT_IF(primaryIterationPacketsMultiplier, "IteratePrimary|includePrimaryEmission")
        ATTRIBUTE_DISPLAYED_IF(primaryIterationPacketsMultiplier, "Level3")

        PROPERTY_DOUBLE(primaryIterationInitialPacketsFraction,
                        "the fraction of primary iteration photon packets to use in the first iteration")
        ATTRIBUTE_MIN_VALUE(primaryIterationInitialPacketsFraction, "]0")
        ATTRIBUTE_MAX_VALUE(primaryIterationInitialPacketsFraction, "1]")
        ATTRIBUTE_DEFAULT_VALUE(primaryIterationInitialPacketsFraction, "1")
        ATTRIBUTE_RELEVANT_IF(primaryIterationInitialPacketsFraction, "IteratePrimary")
        ATTRIBUTE_DISPLAYED_IF(primaryIterationInitialPacketsFraction, "Level3")

        PROPERTY_DOUBLE(primaryIterationPacketsRamp,
                        "the geometric ramp-up factor per primary iteration for the photon packet count")
        ATTRIBUTE_MIN_VALUE(primaryIterationPacketsRamp, "]1")
        ATTRIBUTE_MAX_VALUE(primaryIterationPacketsRamp, "100]")
        ATTRIBUTE_DEFAULT_VALUE(primaryIterationPacketsRamp, "2")
        ATTRIBUTE_RELEVANT_IF(primaryIterationPacketsRamp, "IteratePrimary&primaryIterationInitialPacketsFraction")
        ATTRIBUTE_DISPLAYED_IF(primaryIterationPacketsRamp, "Level3")

        PROPERTY_DOUBLE(primaryIterationPacketsPerCell,
                        "the number of packets per cell for each primary iteration, or zero to disable")
        ATTRIBUTE_MIN_VALUE(primaryIterationPacketsPerCell, "[0")
        ATTRIBUTE_MAX_VALUE(primaryIterationPacketsPerCell, "1e9]")
        ATTRIBUTE_DEFAULT_VALUE(primaryIterationPacketsPerCell, "0")
        ATTRIBUTE_RELEVANT_IF(primaryIterationPacketsPerCell, "IteratePrimary")
        ATTRIBUTE_DISPLAYED_IF(primaryIterationPacketsPerCell, "Level3")

        PROPERTY_DOUBLE(secondaryIterationPacketsMultiplier,
                        "the multiplier on the number of packets for each secondary iteration")
        ATTRIBUTE_MIN_VALUE(secondaryIterationPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(secondaryIterationPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(secondaryIterationPacketsMultiplier, "1")
        ATTRIBUTE_RELEVANT_IF(secondaryIterationPacketsMultiplier, "IterateSecondary")
        ATTRIBUTE_DISPLAYED_IF(secondaryIterationPacketsMultiplier, "Level3")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
