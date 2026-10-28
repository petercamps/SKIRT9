/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TreeSpatialGrid.hpp"
#include "Array.hpp"
#include "ColumnOutFile.hpp"
#include "Configuration.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PathSegmentGenerator.hpp"
#include "Prefetch.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "SpatialGridPlotFile.hpp"
#include "StringUtils.hpp"
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

        // subdivide the nodes that have been flagged, appending their children to the list; the reference to the
        // parent node remains valid while the children are being appended, because appending to a deque never moves
        // its existing elements
        size_t numDivideNodes = divide.sum();
        log->infoSetElapsed(numDivideNodes);
        size_t numDone = 0;
        int numChildren = this->numChildren();
        for (size_t l = 0; l != numEvalNodes; ++l)
        {
            if (divide[l])
            {
                Node& node = nodes[lbeg + l];
                node.setChild(static_cast<int>(nodes.size()));
                for (int c = 0; c != numChildren; ++c)
                    nodes.emplace_back(childExtent(node.extent(), level, c), level + 1);
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
    linkNeighbors();

    // determine the cell indices, the leaf node for each cell, and the initial level of each cell
    int numNodes = static_cast<int>(_nodev.size());
    for (int n = 0; n != numNodes; ++n)
    {
        if (_nodev[n].isLeaf())
        {
            _nodev[n].setCell(static_cast<int>(_idv.size()));
            _idv.push_back(n);
            _initialLevelv.push_back(static_cast<uint8_t>(_nodev[n].level()));
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

void TreeSpatialGrid::linkNeighbors()
{
    // the nodes are in parent-first order, so that the links of each parent are known before its children are linked
    int numNodes = static_cast<int>(_nodev.size());
    for (int p = 0; p != numNodes; ++p)
        if (!_nodev[p].isLeaf()) linkChildren(_nodev, p);
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

namespace
{
    using Node = TreeSpatialGrid::Node;

    // This struct holds the rule for selecting the child of a nonleaf node in a binary tree that contains a given
    // position: the child on the lower or upper side of the splitting plane perpendicular to the axis of the node.
    // A position on the splitting plane is assigned to the upper child.
    struct BinTreeRule
    {
        static int child(const Node& node, double x, double y, double z)
        {
            const double r[3] = {x, y, z};
            int a = node.axis();
            return node.child() + (r[a] < node.center(a) ? 0 : 1);
        }
    };

    // This struct holds the rule for selecting the child of a nonleaf node in an octtree that contains a given
    // position: the octant on the lower or upper side of each of the three planes through the center of the node.
    // A position on one of these planes is assigned to the octant on the upper side of the plane.
    struct OctTreeRule
    {
        static int child(const Node& node, double x, double y, double z)
        {
            return node.child() + (x < node.center(0) ? 0 : 1) + (y < node.center(1) ? 0 : 2)
                   + (z < node.center(2) ? 0 : 4);
        }
    };

    // This function returns the index of the leaf node that contains the specified position, descending the tree from
    // the node with index n in the specified array of nodes, which should contain the position, using the specified
    // rule for selecting a child. If the position is not inside the starting node, the result is not necessarily the
    // leaf node closest to the position, and the caller must verify the result.
    template<class Rule> int descend(const vector<Node>& nodes, int n, double x, double y, double z)
    {
        while (!nodes[n].isLeaf()) n = Rule::child(nodes[n], x, y, z);
        return n;
    }

    // This function returns the index of the leaf node that contains the specified position by searching the tree
    // top-down from its root node, or -1 if the position is outside of the domain.
    template<class Rule> int locate(const vector<Node>& nodes, double x, double y, double z)
    {
        if (!nodes[0].contains(x, y, z)) return -1;
        return descend<Rule>(nodes, 0, x, y, z);
    }

    // This function returns the index of the leaf node that contains the specified position, which must be located
    // just across the specified wall (0-5) of the leaf node with index n. It proceeds from the neighbor across that
    // wall, and returns -1 if the position is outside of the domain. If the neighbor does not contain the position due
    // to numerical inaccuracies (for example, when a path passes very close to the edge or corner of a node), it uses
    // a top-down search as a fall-back.
    template<class Rule> int locateAcross(const vector<Node>& nodes, int n, int wall, double x, double y, double z)
    {
        int m = nodes[n].neighbor(wall);
        if (m >= 0)
        {
            m = descend<Rule>(nodes, m, x, y, z);
            if (nodes[m].contains(x, y, z)) return m;
        }
        return locate<Rule>(nodes, x, y, z);
    }

    // This class implements the path segment generator for a tree grid with the specified rule for selecting a child,
    // as described in the class header.
    template<class Rule> class SegmentGenerator : public PathSegmentGenerator
    {
        Box _extent;                    // the spatial domain of the grid
        double _eps{0.};                // a small distance relative to the extent of the grid
        const vector<Node>& _nodes;     // the array of nodes in the tree
        const vector<int>& _cellNodes;  // the index of the leaf node for each cell
        int _n{-1};                     // index of the leaf node containing the current position

        // the distance to the wall of a node that the path can cross along each axis is
        // (wall - position) * inverse + bias; these values depend only on the direction of the path,
        // and are set up once for each path (see initializeDirection())
        std::array<int, 3> _wallv{{1, 3, 5}};        // the wall (0-5) of a node through which the path can leave it
        std::array<double, 3> _invv{{0., 0., 0.}};   // the reciprocal of the direction component
        std::array<double, 3> _biasv{{0., 0., 0.}};  // zero, or a huge distance for a component that is zero

    public:
        SegmentGenerator(const Box& extent, double eps, const vector<Node>& nodes, const vector<int>& cellNodes)
            : _extent(extent), _eps(eps), _nodes(nodes), _cellNodes(cellNodes)
        {}

        bool next() override
        {
            switch (state())
            {
                case State::KnownCell:
                {
                    // if the initial position is inside the leaf node for the known initial cell, farther from its
                    // walls than a small margin, start from that node without searching, and determine the first
                    // segment
                    int n = _cellNodes[initialCellIndex()];
                    if (_nodes[n].containsWithMargin(rx(), ry(), rz(), _eps))
                    {
                        _n = n;
                        initializeDirection();
                        setState(State::Inside);
                        return nextInside();
                    }
                }

                // otherwise, search for the initial cell as usual
                // intentionally falls through
                case State::Unknown:
                {
                    // try moving the photon packet inside the grid; if this is impossible, return an empty path
                    if (!moveInside(_extent, _eps)) return false;

                    // get the node containing the current location
                    _n = locate<Rule>(_nodes, rx(), ry(), rz());
                    initializeDirection();

                    // if the photon packet started outside the grid, return the corresponding nonzero-length
                    // segment; otherwise fall through to determine the first actual segment
                    if (ds() > 0.) return true;
                }

                // intentionally falls through
                case State::Inside:
                {
                    return nextInside();
                }

                case State::Outside:
                {
                }
            }
            return false;
        }

    private:
        // This function sets up the quantities that depend only on the direction of the path.
        void initializeDirection()
        {
            // The path can leave a node through the upper wall along each axis if the direction component
            // is positive, or through the lower wall otherwise, at a distance (wall - position) / component.
            // To avoid a division and any special cases in each step, we calculate the reciprocal of each
            // component just once for the path. A path parallel to an axis never crosses the walls
            // perpendicular to that axis. We represent this by a zero reciprocal and a huge offset, so that
            // the distance to such a wall is a constant that is never the shortest one. This avoids dividing
            // by zero (which would yield a negative infinity for a negative zero component and not-a-number
            // for a position on the wall).
            const std::array<double, 3> kv{{kx(), ky(), kz()}};
            for (int a = 0; a != 3; ++a)
            {
                bool nonzero = fabs(kv[a]) > 1e-300;
                _wallv[a] = 2 * a + (kv[a] > 0. ? 1 : 0);
                _invv[a] = nonzero ? 1. / kv[a] : 0.;
                _biasv[a] = nonzero ? 0. : DBL_MAX;
            }
        }

        // This function determines the segment from the current position to the first wall of the current node,
        // adjusts the position and the current node accordingly, and returns true.
        bool nextInside()
        {
            const Node& node = _nodes[_n];

            // The next node is one of the three nodes across the walls through which the path can leave
            // this node. Fetching a node from memory takes long compared to the calculations below, and the
            // nodes of a large tree are usually not in the cache. So we ask for all three candidates to be
            // loaded right now, so that this wait overlaps with the calculations. A wall on the boundary of
            // the domain has no neighbor (-1), in which case we harmlessly prefetch the root node instead.
            for (int a = 0; a != 3; ++a) Prefetch::object(&_nodes[std::max(node.neighbor(_wallv[a]), 0)]);

            double dsx = (node.wall(_wallv[0]) - rx()) * _invv[0] + _biasv[0];
            double dsy = (node.wall(_wallv[1]) - ry()) * _invv[1] + _biasv[1];
            double dsz = (node.wall(_wallv[2]) - rz()) * _invv[2] + _biasv[2];

            double ds;
            int wall;
            if (dsx <= dsy && dsx <= dsz)
            {
                ds = dsx;
                wall = _wallv[0];
            }
            else if (dsy <= dsx && dsy <= dsz)
            {
                ds = dsy;
                wall = _wallv[1];
            }
            else
            {
                ds = dsz;
                wall = _wallv[2];
            }
            propagater(ds + _eps);
            setSegment(node.cell(), ds);

            // find the new node by following the link across the crossed wall; this should not fail unless
            // the new location is outside the grid, however on rare occasions it fails due to rounding
            // errors (e.g. in a corner), in which case the function falls back to top-down search
            int oldn = _n;
            _n = locateAcross<Rule>(_nodes, oldn, wall, rx(), ry(), rz());

            // if we're stuck in the same node,
            // try to escape by advancing the position to the next representable coordinates
            if (_n == oldn)
            {
                propagateToNextAfter();
                _n = locate<Rule>(_nodes, rx(), ry(), rz());
            }

            // if we're outside the domain or still stuck in the same node, terminate the path
            if (_n < 0 || _n == oldn) setState(State::Outside);
            return true;
        }
    };
}

////////////////////////////////////////////////////////////////////

int TreeSpatialGrid::cellIndex(Position bfr) const
{
    int n = numChildren() == 2 ? locate<BinTreeRule>(_nodev, bfr.x(), bfr.y(), bfr.z())
                               : locate<OctTreeRule>(_nodev, bfr.x(), bfr.y(), bfr.z());
    return n >= 0 ? _nodev[n].cell() : -1;
}

////////////////////////////////////////////////////////////////////

std::unique_ptr<PathSegmentGenerator> TreeSpatialGrid::createPathSegmentGenerator() const
{
    if (numChildren() == 2) return std::make_unique<SegmentGenerator<BinTreeRule>>(extent(), _eps, _nodev, _idv);
    return std::make_unique<SegmentGenerator<OctTreeRule>>(extent(), _eps, _nodev, _idv);
}

////////////////////////////////////////////////////////////////////

int TreeSpatialGrid::cellLevel(int m) const
{
    return cellNode(m).level();
}

////////////////////////////////////////////////////////////////////

int TreeSpatialGrid::initialCellLevel(int m) const
{
    return _initialLevelv[m];
}

////////////////////////////////////////////////////////////////////

double TreeSpatialGrid::cellSize(int m) const
{
    Box box = cellNode(m).extent();
    return max({box.xwidth(), box.ywidth(), box.zwidth()});
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::faceNeighbors(int m, vector<int>& neighbors) const
{
    neighbors.clear();
    const Node& node = cellNode(m);
    int numChildren = this->numChildren();

    // returns true if the specified node overlaps the wall of the cell perpendicular to axis a in the other two
    // directions, with a nonzero area (so that nodes touching only an edge or a corner of the wall are excluded)
    auto overlaps = [&node](const Node& other, int a) {
        for (int b = 0; b != 3; ++b)
        {
            if (b != a && (other.wall(2 * b) >= node.wall(2 * b + 1) || other.wall(2 * b + 1) <= node.wall(2 * b)))
                return false;
        }
        return true;
    };

    vector<int> stack;
    for (int w = 0; w != 6; ++w)
    {
        int q = node.neighbor(w);
        if (q < 0) continue;

        // the side of the neighbor facing the cell is its lower wall if w is the upper wall of the cell, and vice versa;
        // a child of a node touches that side if its wall on that side coincides with the wall of the node (children
        // copy the walls of their parent, so the coordinates are identical)
        int a = w / 2;
        int facing = w % 2 ? 2 * a : 2 * a + 1;
        stack.push_back(q);
        while (!stack.empty())
        {
            const Node& other = _nodev[stack.back()];
            stack.pop_back();
            if (other.isLeaf())
            {
                if (overlaps(other, a)) neighbors.push_back(other.cell());
            }
            else
            {
                for (int c = 0; c != numChildren; ++c)
                {
                    const Node& child = _nodev[other.child() + c];
                    if (child.wall(facing) == other.wall(facing) && overlaps(child, a))
                        stack.push_back(other.child() + c);
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

vector<int> TreeSpatialGrid::subdivide(const vector<int>& cells)
{
    int numChildren = this->numChildren();
    vector<int> parents;
    parents.reserve(cells.size() * (numChildren - 1));

    for (int m : cells)
    {
        // get the properties of the leaf node holding the cell (by index, because appending invalidates references)
        int p = _idv[m];
        if (!_nodev[p].isLeaf()) throw FATALERROR("Cannot subdivide a cell that is not a leaf node");
        Box extent = _nodev[p].extent();
        int level = _nodev[p].level();

        // the node and cell indices are 32-bit integers
        if (_nodev.size() + numChildren >= static_cast<size_t>(std::numeric_limits<int>::max()))
            throw FATALERROR("The spatial tree grid has too many nodes");

        // append the children and turn the leaf node into a nonleaf node
        int c0 = static_cast<int>(_nodev.size());
        for (int c = 0; c != numChildren; ++c) _nodev.emplace_back(childExtent(extent, level, c), level + 1);
        _nodev[p].setChild(c0);
        _nodev[p].setCell(-1);

        // establish the neighbor links of the children from those of the parent
        linkChildren(_nodev, p);

        // the first child takes over the cell index, and the other children receive new cell indices
        _nodev[c0].setCell(m);
        _idv[m] = c0;
        for (int c = 1; c != numChildren; ++c)
        {
            _nodev[c0 + c].setCell(static_cast<int>(_idv.size()));
            _idv.push_back(c0 + c);
            _initialLevelv.push_back(_initialLevelv[m]);
            parents.push_back(m);
        }
    }
    return parents;
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::reserveCells(int numCells)
{
    size_t numAdded = numCells > this->numCells() ? numCells - this->numCells() : 0;
    size_t numChildren = this->numChildren();
    size_t numSubdivisions = (numAdded + numChildren - 2) / (numChildren - 1);
    _nodev.reserve(_nodev.size() + numSubdivisions * numChildren);
    _idv.reserve(numCells);
    _initialLevelv.reserve(numCells);
}

////////////////////////////////////////////////////////////////////

double TreeSpatialGrid::bytesPerAddedCell() const
{
    double numChildren = this->numChildren();
    return numChildren / (numChildren - 1.) * sizeof(Node) + sizeof(int) + sizeof(uint8_t);
}

////////////////////////////////////////////////////////////////////

void TreeSpatialGrid::writeTopology(ColumnOutFile* outfile) const
{
    outfile->setLongDescription("Topology for tree spatial grid with " + std::to_string(numCells()) + " cells");
    outfile->addColumn("number of children per nonleaf node, followed by the subdivision flag for each node", "", 'd');
    outfile->writeRow(_nodev[0].isLeaf() ? 0 : numChildren());  // zero if the root is not subdivided

    // write the subdivision flags in depth-first order, visiting the children of a node in order of their index;
    // the stack holds the indices of the nodes still to be visited, with the next one to be visited last
    int numChildren = this->numChildren();
    vector<int> stack{0};
    while (!stack.empty())
    {
        const Node& node = _nodev[stack.back()];
        stack.pop_back();
        if (node.isLeaf())
            outfile->writeRow(0);
        else
        {
            outfile->writeRow(1);
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
