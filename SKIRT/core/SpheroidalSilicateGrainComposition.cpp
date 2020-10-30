/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpheroidalSilicateGrainComposition.hpp"
#include "FatalError.hpp"

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

    if (_tableType == TableType::Builtin)
    {
        throw FATALERROR("Spheroidal tables are not part of the SKIRT resources yet!");
        resource = true;
        interpol = _alignmentFraction;
        tableName1 = "NAME OF PRECOMPUTED NON-ALIGNED TABLE";
        tableName2 = "NAME OF PRECOMPUTED ALIGNED TABLE";
    }
    else if (_tableType == TableType::OneTable)
    {
        resource = false;
        interpol = 0;
        tableName1 = _emissionTable;
        tableName2 = string();
    }
    else if (_tableType == TableType::TwoTables)
    {
        resource = false;
        interpol = _alignmentFraction;
        tableName1 = _nonAlignedEmissionTable;
        tableName2 = _alignedEmissionTable;
    }

    return true;
}

////////////////////////////////////////////////////////////////////
