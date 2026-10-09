/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TopologyTreePolicy.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "System.hpp"
#include "TreeNodeEvaluation.hpp"
#include "TreeSpatialGrid.hpp"
#include <fstream>

////////////////////////////////////////////////////////////////////

void TopologyTreePolicy::setupSelfBefore()
{
    TreePolicy::setupSelfBefore();

    // get the information on the grid; don't setup the grid because we are part of it
    _grid = find<TreeSpatialGrid>(false);
    _extent = _grid->extent();

    // open the input file
    string filepath = find<FilePaths>()->input(filename());
    std::ifstream infile = System::ifstream(filepath);
    if (!infile) throw FATALERROR("Could not open the tree topology data file " + filepath);
    auto log = find<Log>();
    log->info(type() + " reads tree topology from text file " + filepath + "...");

    // skip any header lines
    string line;
    while (infile.peek() == '#') getline(infile, line);

    // verify the number of children
    int numChildren = -1;
    infile >> numChildren;
    if (numChildren != 0 && numChildren != 2 && numChildren != 8)
        throw FATALERROR("Tree topology data file has improper format and/or missing data");
    if (numChildren != 0 && numChildren != _grid->numChildren())
        throw FATALERROR("Tree topology data file describes a tree with " + std::to_string(numChildren)
                         + " children per node, while the spatial grid has " + std::to_string(_grid->numChildren())
                         + "; use " + (numChildren == 2 ? "BinTreeSpatialGrid" : "OctTreeSpatialGrid") + " instead");

    // read the subdivision flags in depth-first order into the reference tree; the stack holds the index and level
    // of the nodes for which the flag has not yet been read, with the next node in the file order last
    _childv.assign(1, -1);
    int depth = 0;
    vector<std::pair<int, int>> stack{{0, 0}};
    while (!stack.empty())
    {
        auto [n, level] = stack.back();
        stack.pop_back();
        depth = max(depth, level);

        int flag = -1;
        infile >> flag;
        if (flag == 1 && numChildren != 0)
        {
            int firstChild = static_cast<int>(_childv.size());
            _childv[n] = firstChild;
            _childv.resize(firstChild + numChildren, -1);
            for (int c = numChildren - 1; c >= 0; --c) stack.emplace_back(firstChild + c, level + 1);
        }
        else if (flag != 0)
            throw FATALERROR("Tree topology data file has improper format and/or missing data");
    }
    log->info("Done reading tree topology with " + std::to_string(_childv.size()) + " nodes and depth "
              + std::to_string(depth));

    // warn if the grid cannot reproduce the recorded tree because its maximum level is too small
    if (depth > _grid->maxLevel())
        log->warning("The maximum level of the spatial grid (" + std::to_string(_grid->maxLevel())
                     + ") is smaller than the depth of the recorded tree topology (" + std::to_string(depth) + ")");
}

////////////////////////////////////////////////////////////////////

bool TopologyTreePolicy::needsSubdivide(TreeNodeEvaluation& node) const
{
    // descend the reference tree toward the center of the node, up to the level of the node
    Vec center = node.box().center();
    Box extent = _extent;
    int n = 0;
    for (int level = 0; level != node.level(); ++level)
    {
        // if the reference tree is less deep at this location, the recorded tree does not subdivide the node
        if (_childv[n] < 0) return false;

        int c = _grid->childIndex(extent, level, center);
        extent = _grid->childExtent(extent, level, c);
        n = _childv[n] + c;
    }

    // the recorded tree subdivides the node if the reference node at its level has children
    return _childv[n] >= 0;
}

////////////////////////////////////////////////////////////////////
