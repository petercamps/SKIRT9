/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ResolvedSpheresTreePolicy.hpp"
#include "Array.hpp"
#include "Log.hpp"
#include "TextInFile.hpp"
#include "TreeNodeEvaluation.hpp"

////////////////////////////////////////////////////////////////////

void ResolvedSpheresTreePolicy::setupSelfBefore()
{
    TreePolicy::setupSelfBefore();

    // read the spheres, applying the configured number of bins and reach unless overridden
    TextInFile infile(this, filename(), "spheres");
    infile.addColumn("position x", "length", "pc");
    infile.addColumn("position y", "length", "pc");
    infile.addColumn("position z", "length", "pc");
    infile.addColumn("radius", "length", "pc");
    if (importNumBins()) infile.addColumn("number of bins");
    if (importReach()) infile.addColumn("reach");
    Array row;
    while (infile.readRow(row))
    {
        double radius = row[3];
        if (radius > 0.)
        {
            int i = 4;
            double numBins = importNumBins() && row[i] > 0. ? row[i] : this->numBins();
            if (importNumBins()) i++;
            double reach = importReach() && row[i] > 0. ? row[i] : this->reach();

            _centerv.emplace_back(row[0], row[1], row[2]);
            _reachv.push_back(reach * radius);
            _sizev.push_back(radius / numBins);
        }
    }
    infile.close();
    int numSpheres = _centerv.size();
    find<Log>()->info(type() + " uses " + std::to_string(numSpheres) + " spheres");

    // organize the enlarged spheres in a search structure
    auto bounds = [this](int m) {
        Vec r(_reachv[m], _reachv[m], _reachv[m]);
        return Box(_centerv[m] - r, _centerv[m] + r);
    };
    auto intersects = [this](int m, const Box& box) { return box.intersects(_centerv[m], _reachv[m]); };
    _search.loadEntities(numSpheres, bounds, intersects);
}

////////////////////////////////////////////////////////////////////

bool ResolvedSpheresTreePolicy::needsSubdivide(TreeNodeEvaluation& node) const
{
    const Box& box = node.box();
    double width = max({box.xwidth(), box.ywidth(), box.zwidth()});
    for (int m : _search.entitiesFor(box))
    {
        if (width > _sizev[m] && box.intersects(_centerv[m], _reachv[m])) return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////
