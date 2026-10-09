/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TreeSpatialGrid.hpp"
#include "Array.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "SpatialGridPlotFile.hpp"
#include "StringUtils.hpp"
#include "TextOutFile.hpp"
#include "TreeNodeEvaluation.hpp"

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::setupSelfBefore()
{
    BoxSpatialGrid::setupSelfBefore();

    if (maxLevel() < minLevel()) throw FATALERROR("Maximum tree level cannot be below minimum level");
}

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of nodes evaluated between two invocations of infoIfElapsed()
    const size_t logEvalChunkSize = 10000;

    // maximum number of nodes subdivided between two invocations of infoIfElapsed()
    const size_t logDivideChunkSize = 5000;
}

////////////////////////////////////////////////////////////////////

std::deque<TreeSpatialGrid::Node> TreeSpatialGrid::constructTree() const
{
    auto log = find<Log>();
    auto parallel = find<ParallelFactory>()->parallelDistributed();

    // get the information needed for evaluating the properties of the media in a node
    auto ms = find<MediumSystem>(false);  // don't setup the medium system because we are part of it
    auto random = find<Random>();
    int numSamples = find<Configuration>()->numDensitySamples();

    // initialize the node list with the root node; the list is a deque so that adding nodes never moves existing ones
    std::deque<Node> nodes;
    nodes.emplace_back(extent(), 0);

    // initialize iteration variables to level 0
    int level = 0;    // current level
    size_t lbeg = 0;  // node index range for the current level;
    size_t lend = 1;  // at level 0, the node list contains just the root node

    // subdivide nodes level by level until all nodes satisfy the configured criteria
    while (lend != lbeg)
    {
        size_t numEvalNodes = lend - lbeg;
        log->info("Subdividing level " + std::to_string(level) + ": " + std::to_string(numEvalNodes) + " nodes");
        log->infoSetElapsed(numEvalNodes);

        // evaluate nodes at this level: value in the array becomes one for nodes that need to be subdivided;
        // we parallelize this operation because it might be resource intensive (e.g. sampling densities)
        Array divide(numEvalNodes);
        if (level < minLevel())
        {
            divide = 1.;
        }
        else if (level < maxLevel() && !_policies.empty())
        {
            parallel->call(numEvalNodes, [this, log, ms, random, numSamples, level, lbeg, &nodes,
                                          &divide](size_t firstIndex, size_t numIndices) {
                // each thread uses its own instance for evaluating nodes, reset for each node
                TreeNodeEvaluation evaluation(ms, random, numSamples);
                while (numIndices)
                {
                    size_t currentChunkSize = min(logEvalChunkSize, numIndices);
                    for (size_t l = firstIndex; l != firstIndex + currentChunkSize; ++l)
                    {
                        evaluation.reset(nodes[lbeg + l].extent(), level);
                        for (auto policy : _policies)
                        {
                            if (policy->needsSubdivide(evaluation))
                            {
                                divide[l] = 1.;
                                break;
                            }
                        }
                    }
                    log->infoIfElapsed("Evaluation for level " + std::to_string(level) + ": ", currentChunkSize);
                    firstIndex += currentChunkSize;
                    numIndices -= currentChunkSize;
                }
            });
            ProcessManager::sumToAll(divide);
        }

        // subdivide the nodes that have been flagged, appending their children to the list
        size_t numDivideNodes = divide.sum();
        log->infoSetElapsed(numDivideNodes);
        size_t numDone = 0;
        for (size_t l = 0; l != numEvalNodes; ++l)
        {
            if (divide[l])
            {
                Node& node = nodes[lbeg + l];
                node.setChild(static_cast<int>(nodes.size()));
                appendChildren(node, nodes);
                numDone++;
                if (numDone % logDivideChunkSize == 0)
                    log->infoIfElapsed("Subdivision for level " + std::to_string(level) + ": ", logDivideChunkSize);
            }
        }

        // the node and cell indices are 32-bit integers
        if (nodes.size() >= static_cast<size_t>(std::numeric_limits<int>::max()))
            throw FATALERROR("The spatial tree grid has too many nodes");

        // update iteration variables to the next level
        level++;
        lbeg = lend;
        lend = nodes.size();
    }
    return nodes;
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::setupSelfAfter()
{
    BoxSpatialGrid::setupSelfAfter();

    // determine a small fraction relative to the spatial extent of the grid; used during path traversal
    _eps = 1e-12 * extent().diagonal();

    // construct the tree and copy it into a contiguous array, releasing the memory held by the construction list
    Log* log = find<Log>();
    log->info("Constructing the spatial tree grid...");
    {
        std::deque<Node> nodes = constructTree();
        _nodev.assign(nodes.cbegin(), nodes.cend());
    }

    // establish the neighbor links
    linkNeighbors(_nodev);

    // determine the cell indices and the leaf node for each cell
    int numNodes = static_cast<int>(_nodev.size());
    for (int n = 0; n != numNodes; ++n)
    {
        if (_nodev[n].isLeaf())
        {
            _nodev[n].setCell(static_cast<int>(_idv.size()));
            _idv.push_back(n);
        }
    }

    // determine the number of cells at each level in the tree hierarchy
    vector<int> countv;
    int numCells = _idv.size();
    for (int m = 0; m != numCells; ++m)
    {
        int level = cellNode(m).level();
        if (level + 1 > static_cast<int>(countv.size())) countv.resize(level + 1);
        countv[level]++;
    }

    // log these statistics, including a basic histogram
    log->info("Finished construction of the spatial tree grid");
    log->info("Number of cells at each level in the tree hierarchy:");
    int numLevels = countv.size();
    int maxCount = *std::max_element(countv.cbegin(), countv.cend());
    for (int level = 0; level != numLevels; ++level)
    {
        size_t numStars = std::round(20. * countv[level] / maxCount);
        log->info("  Level " + StringUtils::toString(level, 'd', 0, 2) + ":"
                  + StringUtils::toString(countv[level], 'd', 0, 9) + " ("
                  + StringUtils::toString(100. * countv[level] / numCells, 'f', 1, 5) + "%)  |"
                  + string(numStars, '*'));
    }
    log->info("  TOTAL   :" + StringUtils::toString(numCells, 'd', 0, 9) + " (100.0%)");

    // make the BoxCellDensityMixIn verify whether we can offer the DensityInCellInterface
    BoxCellDensityMixIn::setup(this);
}

////////////////////////////////////////////////////////////////////

int TreeSpatialGrid::numCells() const
{
    return _idv.size();
}

////////////////////////////////////////////////////////////////////

Box TreeSpatialGrid::cellBox(int m) const
{
    return cellNode(m).extent();
}

////////////////////////////////////////////////////////////////////

double TreeSpatialGrid::volume(int m) const
{
    return cellNode(m).extent().volume();
}

////////////////////////////////////////////////////////////////////

double TreeSpatialGrid::diagonal(int m) const
{
    return cellNode(m).extent().diagonal();
}

////////////////////////////////////////////////////////////////////

Position TreeSpatialGrid::centralPositionInCell(int m) const
{
    return Position(cellNode(m).extent().center());
}

////////////////////////////////////////////////////////////////////

Position TreeSpatialGrid::randomPositionInCell(int m) const
{
    return random()->position(cellNode(m).extent());
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::writeTopology(TextOutFile* outfile) const
{
    outfile->writeLine("# Topology for tree spatial grid with " + std::to_string(numCells()) + " cells");
    outfile->writeLine(std::to_string(_nodev[0].isLeaf() ? 0 : numChildren()));  // zero if the root is not subdivided

    // write the subdivision flags in depth-first order, visiting the children of a node in order of their index;
    // the stack holds the indices of the nodes still to be visited, with the next one to be visited last
    int numChildren = this->numChildren();
    vector<int> stack{0};
    while (!stack.empty())
    {
        const Node& node = _nodev[stack.back()];
        stack.pop_back();
        if (node.isLeaf())
            outfile->writeLine("0");
        else
        {
            outfile->writeLine("1");
            for (int c = numChildren - 1; c >= 0; --c) stack.push_back(node.child() + c);
        }
    }
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    // Output the root cell and all leaf cells that are close to the section plane
    outfile->writeRectangle(xmin(), ymin(), xmax(), ymax());
    int nCells = numCells();
    for (int m = 0; m != nCells; ++m)
    {
        Box box = cellBox(m);
        if (fabs(box.zmin()) < 1e-8 * extent().zwidth())
        {
            outfile->writeRectangle(box.xmin(), box.ymin(), box.xmax(), box.ymax());
        }
    }
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
{
    // Output the root cell and all leaf cells that are close to the section plane
    outfile->writeRectangle(xmin(), zmin(), xmax(), zmax());
    int nCells = numCells();
    for (int m = 0; m != nCells; ++m)
    {
        Box box = cellBox(m);
        if (fabs(box.ymin()) < 1e-8 * extent().ywidth())
        {
            outfile->writeRectangle(box.xmin(), box.zmin(), box.xmax(), box.zmax());
        }
    }
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::write_yz(SpatialGridPlotFile* outfile) const
{
    // Output the root cell and all leaf cells that are close to the section plane
    outfile->writeRectangle(ymin(), zmin(), ymax(), zmax());
    int nCells = numCells();
    for (int m = 0; m != nCells; ++m)
    {
        Box box = cellBox(m);
        if (fabs(box.xmin()) < 1e-8 * extent().xwidth())
        {
            outfile->writeRectangle(box.ymin(), box.zmin(), box.ymax(), box.zmax());
        }
    }
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::write_xyz(SpatialGridPlotFile* outfile) const
{
    // determine the number of cells at each level in the tree hierarchy
    vector<int> countv;
    int nCells = numCells();
    for (int m = 0; m != nCells; ++m)
    {
        int level = cellNode(m).level();
        if (level + 1 > static_cast<int>(countv.size())) countv.resize(level + 1);
        countv[level]++;
    }

    // determine the number of levels to be included in output
    int maxLevel = static_cast<int>(countv.size()) - 1;
    int highestWriteLevel = 0;
    int cumulativeCells = 0;
    for (; highestWriteLevel <= maxLevel; ++highestWriteLevel)
    {
        cumulativeCells += countv[highestWriteLevel];
        if (cumulativeCells > 2500) break;  // empirical number
    }

    // inform the user if we are limiting output
    if (highestWriteLevel < maxLevel)
        find<Log>()->info("Limiting 3D grid plot output tree to level " + std::to_string(highestWriteLevel) + ", i.e. "
                          + std::to_string(cumulativeCells) + " cells.");

    // output all leaf cells up to a certain level
    for (int m = 0; m != nCells; ++m)
    {
        const Node& node = cellNode(m);
        if (node.level() <= highestWriteLevel)
        {
            Box box = node.extent();
            outfile->writeCube(box.xmin(), box.ymin(), box.zmin(), box.xmax(), box.ymax(), box.zmax());
        }
    }
}

////////////////////////////////////////////////////////////////////

bool TreeSpatialGrid::offersInterface(const std::type_info& interfaceTypeInfo) const
{
    if (interfaceTypeInfo == typeid(DensityInCellInterface)) return BoxCellDensityMixIn::offersInterface();
    return BoxSpatialGrid::offersInterface(interfaceTypeInfo);
}

////////////////////////////////////////////////////////////////////
