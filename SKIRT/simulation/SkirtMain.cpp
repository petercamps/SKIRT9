/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BuildInfo.hpp"
#include "ProcessManager.hpp"
#include "SignalHandler.hpp"
#include "SimulationItemRegistry.hpp"
#include "SkirtCommandLineHandler.hpp"
#include "System.hpp"

//////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    // Initialize inter-process communication capability, if present
    ProcessManager pm(&argc, &argv);

    // Initialize the system and install signal handlers
    System system(argc, argv);
    SignalHandler::InstallSignalHandlers();

    // Add all simulation items to the item registry. The second argument is the version of the ski file format, which
    // changes only with an incompatible change to the structure of the file, and not with changes to the simulation
    // items and their properties, which the ski file upgrade procedure of PTS handles. The format has been "9" since
    // SKIRT 9, and does not follow the SKIRT version.
    string version = BuildInfo::projectVersion();
    SimulationItemRegistry registry(version, "9");

    // handle the command line arguments
    SkirtCommandLineHandler handler;
    return handler.perform();
}

//////////////////////////////////////////////////////////////////////
