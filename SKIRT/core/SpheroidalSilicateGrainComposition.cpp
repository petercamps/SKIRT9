/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpheroidalSilicateGrainComposition.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

void SpheroidalSilicateGrainComposition::setupSelfBefore()
{

    PolarizedSilicateGrainComposition::setupSelfBefore();

    if (!_alignedEmissionTable.empty() && _alignmentFraction != 1. && _nonAlignedEmissionTable.empty())
    {
        throw FATALERROR(
            "No custom emission table provided for non-aligned grains, yet the alignment fraction is not 1!");
    }
}

////////////////////////////////////////////////////////////////////

string SpheroidalSilicateGrainComposition::name() const
{
    return "Spheroidal_Polarized_Draine_Silicate";
}

////////////////////////////////////////////////////////////////////

bool SpheroidalSilicateGrainComposition::resourcesForSpheroidalEmission(bool& resource, double& interpol,
                                                                        std::string& tableName1,
                                                                        std::string& tableName2) const
{

    if (!_alignedEmissionTable.empty())
    {
        resource = false;
        if (!_nonAlignedEmissionTable.empty())
        {
            interpol = _alignmentFraction;
            tableName1 = _nonAlignedEmissionTable;
            tableName2 = _alignedEmissionTable;
        }
        else
        {
            interpol = 0;
            tableName1 = _alignedEmissionTable;
            tableName2 = string();
        }
    }
    else
    {
        throw FATALERROR("Spheroidal tables are not part of the SKIRT resources yet!");
        resource = true;
        interpol = _alignmentFraction;
        tableName1 = "NAME OF PRECOMPUTED NON-ALIGNED TABLE";
        tableName2 = "NAME OF PRECOMPUTED ALIGNED TABLE";
    }

    return true;
}

////////////////////////////////////////////////////////////////////
