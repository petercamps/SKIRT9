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

    ENUM_DEF(TableType, Builtin, OneTable, TwoTables)
        ENUM_VAL(TableType, Builtin, "builtin resources")
        ENUM_VAL(TableType, OneTable, "single custom table")
        ENUM_VAL(TableType, TwoTables, "two custom tables with interpolation")
    ENUM_END()

    ITEM_CONCRETE(SpheroidalSilicateGrainComposition, PolarizedSilicateGrainComposition,
                  "a spheroidal silicate dust grain composition with support for polarization")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SpheroidalSilicateGrainComposition, "Spheroidal")

        PROPERTY_ENUM(tableType, TableType, "the type of emission tables to use")
        ATTRIBUTE_DEFAULT_VALUE(tableType, "Builtin")

        PROPERTY_STRING(emissionTable, "the name of the file tabulating properties for polarized emission by "
                                       "arbitrarily aligned spheroidal grains")
        ATTRIBUTE_RELEVANT_IF(emissionTable, "tableTypeOneTable")

        PROPERTY_STRING(
            alignedEmissionTable,
            "the name of the file tabulating properties for polarized emission by perfectly aligned spheroidal grains")
        ATTRIBUTE_RELEVANT_IF(alignedEmissionTable, "tableTypeTwoTables")

        PROPERTY_STRING(
            nonAlignedEmissionTable,
            "the name of the file tabulating properties for polarized emission by non-aligned spheroidal grains")
        ATTRIBUTE_RELEVANT_IF(nonAlignedEmissionTable, "tableTypeTwoTables")

        PROPERTY_DOUBLE(alignmentFraction,
                        "the alignment fraction of the spheroidal grains with the local magnetic field")
        ATTRIBUTE_DEFAULT_VALUE(alignmentFraction, "1.")
        ATTRIBUTE_MIN_VALUE(alignmentFraction, "0.")
        ATTRIBUTE_MAX_VALUE(alignmentFraction, "1.")
        ATTRIBUTE_RELEVANT_IF(alignmentFraction, "tableTypeBuiltin|tableTypeTwoTables")

    ITEM_END()

public:
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
