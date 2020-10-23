/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHEROIDALSILICATEGRAINCOMPOSITION_HPP
#define SPHEROIDALSILICATEGRAINCOMPOSITION_HPP

#include "PolarizedSilicateGrainComposition.hpp"

////////////////////////////////////////////////////////////////////

/** The SpheroidalSilicateGrainComposition class represents the optical and calorimetric properties
    of spheroidal silicate dust grains with partial support for polarization. More precisely, the
    current implementation supports polarized thermal emission by (partially) aligned spheroidal
    grains, but assumes spherical grains for scattering and absorption interactions.

    The optical scattering and absorption properties and the calorimetric properties are taken from
    the PolarizedSilicateGrainComposition class, from which this class derives. The optical
    properties driving the polarization signature for therrmal emission are obtained from
    additional built-in tables or can be provided by the user, as described below.

    TO DO: describe built-in and file input options. */
class SpheroidalSilicateGrainComposition : public PolarizedSilicateGrainComposition
{
    ITEM_CONCRETE(SpheroidalSilicateGrainComposition, PolarizedSilicateGrainComposition,
                  "a spheroidal silicate dust grain composition with support for polarization")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SpheroidalSilicateGrainComposition, "Spheroidal")

        PROPERTY_STRING(
            alignedEmissionTable,
            "the name of the file tabulating properties for polarized emission by perfectly aligned spheroidal grains")
        ATTRIBUTE_DEFAULT_VALUE(alignedEmissionTable, "")
        ATTRIBUTE_REQUIRED_IF(alignedEmissionTable, "false")

        PROPERTY_STRING(
            nonAlignedEmissionTable,
            "the name of the file tabulating properties for polarized emission by non-aligned spheroidal grains")
        ATTRIBUTE_DEFAULT_VALUE(nonAlignedEmissionTable, "")
        ATTRIBUTE_REQUIRED_IF(nonAlignedEmissionTable, "false")

        PROPERTY_DOUBLE(alignmentFraction,
                        "the alignment fraction of the spheroidal grains with the local magnetic field")
        ATTRIBUTE_DEFAULT_VALUE(alignmentFraction, "1.")
        ATTRIBUTE_MIN_VALUE(alignmentFraction, "0.")
        ATTRIBUTE_MAX_VALUE(alignmentFraction, "1.")

    ITEM_END()

public:
    /** This function verifies that all attribute values have been appropriately set. */
    void setupSelfBefore() override;

    /** This function returns a brief human-readable identifier for this grain composition. */
    string name() const override;

    /** This function returns information on the resources required for implementing thermal
        emission from aligned spheriodal grains. For more information, see
        GrainComposition::resourcesForSpheroidalEmission. */
    bool resourcesForSpheroidalEmission(bool& resource, double& interpol, string& tableName1,
                                        string& tableName2) const override;
};

////////////////////////////////////////////////////////////////////

#endif
