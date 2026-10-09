/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SiteListTreePolicy.hpp"
#include "Array.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "SiteListInterface.hpp"
#include "TextInFile.hpp"
#include "TreeNodeEvaluation.hpp"
#include "TreeSpatialGrid.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // a node of the reference tree at the level being constructed, with the range of its sites in the index list
    struct RangeNode
    {
        int node;      // index of the node in the reference tree
        Box extent;    // spatial extent of the node
        size_t begin;  // index in the index list of the first site in the node
        size_t end;    // index in the index list beyond the last site in the node
    };
}

////////////////////////////////////////////////////////////////////

void SiteListTreePolicy::setupSelfBefore()
{
    TreePolicy::setupSelfBefore();

    // get the information on the grid; don't setup the grid because we are part of it
    _grid = find<TreeSpatialGrid>(false);
    _extent = _grid->extent();
    _minLevel = _grid->minLevel();
    int maxLevel = _grid->maxLevel();
    int numChildren = _grid->numChildren();

    // get the site positions, either from the input file or from the first medium component offering a site list
    vector<Vec> sitev;
    if (!filename().empty())
    {
        TextInFile infile(this, filename(), "site positions");
        infile.addColumn("position x", "length", "pc");
        infile.addColumn("position y", "length", "pc");
        infile.addColumn("position z", "length", "pc");
        Array row;
        while (infile.readRow(row)) sitev.emplace_back(row[0], row[1], row[2]);
    }
    else
    {
        // search the media and their geometries; don't setup the medium system because we are part of it
        auto ms = find<MediumSystem>(false);
        auto sli = ms ? ms->interface<SiteListInterface>(0, 2, false) : nullptr;
        if (!sli) throw FATALERROR("There is no input file with site positions and no medium offering a site list");
        int numSites = sli->numSites();
        sitev.reserve(numSites);
        for (int m = 0; m != numSites; ++m) sitev.push_back(sli->sitePosition(m));
    }

    // ignore sites outside of the spatial domain
    sitev.erase(std::remove_if(sitev.begin(), sitev.end(), [this](Vec r) { return !_extent.contains(r); }),
                sitev.end());
    find<Log>()->info(type() + " uses " + std::to_string(sitev.size()) + " sites inside the spatial domain");

    // construct the reference tree level by level, subdividing the nodes with more than one site below the maximum
    // level; the sites in each node at the current level form a consecutive range in a list of site indices, which
    // is reordered when the node is subdivided so that the sites in each child again form a consecutive range
    vector<int> idv(sitev.size());
    std::iota(idv.begin(), idv.end(), 0);
    vector<int> childIndexv;  // the index of the child containing each site in the node being subdivided
    vector<int> sortedv;      // the site indices in the node being subdivided, sorted by child index
    vector<size_t> startv;    // the index in the sorted list of the first site in each child
    vector<size_t> nextv;     // the index in the sorted list of the next site to be stored for each child

    _childv.assign(1, -1);
    vector<RangeNode> levelNodes{{0, _extent, 0, idv.size()}};
    vector<RangeNode> nextNodes;
    for (int level = 0; level != maxLevel && !levelNodes.empty(); ++level)
    {
        nextNodes.clear();
        for (const RangeNode& parent : levelNodes)
        {
            size_t numSites = parent.end - parent.begin;
            if (numSites < 2) continue;

            // subdivide the node
            int firstChild = static_cast<int>(_childv.size());
            _childv[parent.node] = firstChild;
            _childv.resize(firstChild + numChildren, -1);

            // sort the sites in the node by the index of the child containing them (a counting sort)
            childIndexv.resize(numSites);
            startv.assign(numChildren + 1, 0);
            for (size_t i = 0; i != numSites; ++i)
            {
                int c = _grid->childIndex(parent.extent, level, sitev[idv[parent.begin + i]]);
                childIndexv[i] = c;
                startv[c + 1]++;
            }
            for (int c = 0; c != numChildren; ++c) startv[c + 1] += startv[c];
            sortedv.resize(numSites);
            nextv.assign(startv.cbegin(), startv.cend() - 1);
            for (size_t i = 0; i != numSites; ++i) sortedv[nextv[childIndexv[i]]++] = idv[parent.begin + i];
            std::copy(sortedv.cbegin(), sortedv.cend(), idv.begin() + parent.begin);

            // remember the children that again contain more than one site
            for (int c = 0; c != numChildren; ++c)
            {
                if (startv[c + 1] - startv[c] > 1)
                    nextNodes.push_back({firstChild + c, _grid->childExtent(parent.extent, level, c),
                                         parent.begin + startv[c], parent.begin + startv[c + 1]});
            }
        }
        std::swap(levelNodes, nextNodes);
    }
}

////////////////////////////////////////////////////////////////////

bool SiteListTreePolicy::needsSubdivide(TreeNodeEvaluation& node) const
{
    // descend the reference tree toward the center of the node, up to the level of the node or to a leaf
    Vec center = node.box().center();
    Box extent = _extent;
    int n = 0;
    int level = 0;
    while (level != node.level() && _childv[n] >= 0)
    {
        int c = _grid->childIndex(extent, level, center);
        extent = _grid->childExtent(extent, level, c);
        n = _childv[n] + c;
        level++;
    }

    // if the reference tree is still subdivided at the level of the node, the node contains multiple sites
    if (_childv[n] >= 0) return true;

    // otherwise, subdivide up to the extra levels below the leaf of the reference tree; if that leaf is coarser
    // than the minimum level of the grid, which subdivides down to that level anyway, count from the minimum level
    return node.level() < max(level, _minLevel) + numExtraLevels();
}

////////////////////////////////////////////////////////////////////
