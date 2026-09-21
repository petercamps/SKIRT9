/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BinTreeSpatialGrid.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "PathSegmentGenerator.hpp"
#include "Random.hpp"
#include "SpatialGridPlotFile.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <array>
#include <fstream>
#include <limits>

////////////////////////////////////////////////////////////////////

namespace
{
    // The Node class represents a node in the binary tree, which is either a leaf node (corresponding to a spatial
    // cell) or a nonleaf node with two children. The nodes are stored by value in a single contiguous array, and they
    // refer to each other by index into this array rather than by pointer. The two children of a nonleaf node are
    // always consecutive in the array, so a node needs to store just the index of its first child. A node does not
    // have a parent link because it is never needed. The node class has no virtual functions.
    //
    // The six walls of a node are numbered such that wall 2*a is the lower and wall 2*a+1 the upper wall along the
    // axis a, with x=0, y=1 and z=2. In other words, the walls are numbered as BACK, FRONT, LEFT, RIGHT, BOTTOM, TOP
    // in the TreeNode class. A node stores the coordinates of its walls in this order, as well as the index of its
    // neighbor across each of its walls (see the Tree class for details on how the neighbors are defined).
    //
    // A nonleaf node at level l is split at its geometrical center along the axis l%3, so that the coordinates of the
    // resulting children are identical to the ones obtained for the corresponding TreeNode subclass. The lower child
    // (the one with the lowest coordinates along the split axis) comes first.
    //
    // The size of a node is 88 bytes: 6 doubles for the walls and 10 integers for the neighbors, children, cell
    // index, level and axis. The 32-bit integer indices keep the node small.
    class Node
    {
    public:
        // this constructor creates a node with the specified extent and level, without any children or neighbors
        Node(const Box& extent, int level)
            : _wallv{{extent.xmin(), extent.xmax(), extent.ymin(), extent.ymax(), extent.zmin(), extent.zmax()}},
              _level(level), _axis(level % 3)
        {}

        // this function returns the extent of the node as a box
        Box extent() const { return Box(_wallv[0], _wallv[2], _wallv[4], _wallv[1], _wallv[3], _wallv[5]); }

        // this function returns the level of the node in the tree, with level zero for the root node
        int level() const { return _level; }

        // this function returns the axis (0=x, 1=y, 2=z) perpendicular to the plane along which a nonleaf node is split
        int axis() const { return _axis; }

        // this function returns the coordinate of the plane along which a nonleaf node is split (along its axis)
        double split() const { return 0.5 * (_wallv[2 * _axis] + _wallv[2 * _axis + 1]); }

        // this function returns the coordinate of the specified wall (0-5) along the axis perpendicular to that wall
        double wall(int w) const { return _wallv[w]; }

        // this function returns true if the specified position is inside the node, borders included
        bool contains(double x, double y, double z) const
        {
            return x >= _wallv[0] && x <= _wallv[1] && y >= _wallv[2] && y <= _wallv[3] && z >= _wallv[4]
                   && z <= _wallv[5];
        }

        // this function returns true if the node is a leaf node, i.e. a spatial cell without children
        bool isLeaf() const { return _child < 0; }

        // this function returns the index of the first child of a nonleaf node; the second child follows it directly
        int child() const { return _child; }

        // this function returns the cell index for a leaf node, or -1 for a nonleaf node
        int cell() const { return _cell; }

        // this function returns the index of the neighbor across the specified wall (0-5), or -1 if there is none
        int neighbor(int w) const { return _neighborv[w]; }

        // this function returns a new node representing the lower half of this node, without children or neighbors
        Node lowerChild() const
        {
            Node node(extent(), _level + 1);
            node._wallv[2 * _axis + 1] = split();
            return node;
        }

        // this function returns a new node representing the upper half of this node, without children or neighbors
        Node upperChild() const
        {
            Node node(extent(), _level + 1);
            node._wallv[2 * _axis] = split();
            return node;
        }

        // these functions set the index of the first child, the cell index, and the index of the neighbor at a wall
        void setChild(int child) { _child = child; }
        void setCell(int cell) { _cell = cell; }
        void setNeighbor(int w, int neighbor) { _neighborv[w] = neighbor; }

    private:
        std::array<double, 6> _wallv;                             // xmin, xmax, ymin, ymax, zmin, zmax
        std::array<int, 6> _neighborv{{-1, -1, -1, -1, -1, -1}};  // neighbor across each wall
        int _child{-1};                                           // index of the first child; -1 for a leaf node
        int _cell{-1};                                            // cell index; -1 for a nonleaf node
        int _level{0};                                            // level of the node in the tree
        int _axis{0};                                             // axis along which the node is split (level % 3)
    };

    // This function asks the processor to load the memory at the specified address into its cache in the background.
    // It is just a hint without any other effect: it does not wait for the memory to arrive and it can never fail,
    // even for an invalid address. Prefetching is not part of standard C++. The builtin function offering it is
    // available in GCC, Clang and the Intel compilers; for other compilers the hint is skipped.
    // NOTE: this is a deliberate, local exception to the coding style that avoids conditional compilation.
#if defined(__GNUC__) || defined(__clang__)
    inline void prefetch(const void* address)
    {
        __builtin_prefetch(address);
    }
#else
    inline void prefetch(const void*) {}
#endif

    // This function asks the processor to load all cache lines holding the specified node into its cache in the
    // background. A node is larger than 64 bytes, so it can span three cache lines of 64 bytes (and two of 128 bytes).
    // Requesting its first byte, its last byte, and the byte 64 positions after the first one covers all of them.
    static_assert(sizeof(Node) > 64 && sizeof(Node) <= 128, "prefetchNode() must be adjusted to the size of a node");
    inline void prefetchNode(const Node* node)
    {
        const char* first = static_cast<const char*>(static_cast<const void*>(node));
        prefetch(first);
        prefetch(first + 64);
        prefetch(first + sizeof(Node) - 1);
    }

    // This function returns the index of the leaf node that contains the specified position. The search starts at the
    // node with the specified index (which should contain the position) in the specified array of nodes, and descends
    // the tree by selecting the child on the appropriate side of the splitting plane until it reaches a leaf node.
    // This is the same selection rule as in the BinTreeNode class. If the position is not inside the starting node, the
    // result is not necessarily the leaf node closest to the position, and the caller must verify the result.
    int descend(const Node* nodes, int n, double x, double y, double z)
    {
        const double r[3] = {x, y, z};
        const Node* node = nodes + n;
        while (!node->isLeaf())
        {
            n = node->child() + (r[node->axis()] < node->split() ? 0 : 1);
            node = nodes + n;
        }
        return n;
    }

    // This function returns the index of the leaf node that contains the specified position by searching the tree
    // top-down from its root node (the first node in the specified array). It returns -1 if the position is outside of
    // the domain.
    int locate(const Node* nodes, double x, double y, double z)
    {
        if (!nodes[0].contains(x, y, z)) return -1;
        return descend(nodes, 0, x, y, z);
    }

    // This function returns the index of the leaf node that contains the specified position, which must be located
    // just across the specified wall (0-5) of the leaf node with index n. It proceeds from the neighbor across that
    // wall, and returns -1 if the position is outside of the domain. In rare cases, the neighbor does not contain the
    // position due to numerical inaccuracies (for example, when a path passes very close to the edge or corner of a
    // node), so that a top-down search is used as a fall-back.
    int locateAcross(const Node* nodes, int n, int wall, double x, double y, double z)
    {
        int m = nodes[n].neighbor(wall);
        if (m >= 0)
        {
            m = descend(nodes, m, x, y, z);
            if (nodes[m].contains(x, y, z)) return m;
        }
        return locate(nodes, x, y, z);
    }

    // This function reads the tree topology from the specified input stream, and returns the corresponding subdivision
    // flags: 1 if the node is subdivided and 0 if not, in the order of a depth-first traversal of the tree (i.e. the
    // order in which they occur in the file). The function stops reading as soon as the tree is complete.
    vector<char> readTopology(std::ifstream& infile)
    {
        // skip any header lines
        string line;
        while (infile.peek() == '#') getline(infile, line);

        // verify that the number of children per node is acceptable for a binary tree
        int numChildren = -1;
        infile >> numChildren;
        if (!infile) throw FATALERROR("Topology input file has improper format and/or missing data");
        switch (numChildren)
        {
            case 0:
            case 2: break;
            case 8: throw FATALERROR("Topology input file specifies an octtree, which is not supported by this grid");
            default:
                throw FATALERROR("Topology input file specifies unsupported number of children in node subdivision");
        }

        // read a flag for each node, keeping track of the number of nodes still waiting for their flag:
        // each node consumes a flag, and each subdivided node adds two new nodes, for a net gain of one
        vector<char> flagv;
        size_t numPending = 1;
        while (numPending)
        {
            int flag = -1;
            infile >> flag;
            if (!infile || (flag != 0 && flag != 1))
                throw FATALERROR("Topology input file has improper format and/or missing data");
            flagv.push_back(static_cast<char>(flag));
            if (flag)
                numPending++;
            else
                numPending--;
        }

        // the root node can't be subdivided if there are no children
        if (numChildren == 0 && flagv.size() != 1)
            throw FATALERROR("Topology input file specifies subdivision without specifying the number of children");
        return flagv;
    }
}

////////////////////////////////////////////////////////////////////

// The Tree class holds the complete tree data structure, consisting of the array of nodes (leaf and nonleaf)
// and the array that maps cell indices to nodes. The tree is immutable after construction.
//
// The nodes are numbered in the same order as in the other tree grids, so that the cell indices are the same for
// a given topology. That order results from subdividing the tree depth-first: the root node has index zero, and
// each time a node is subdivided, its two children are appended to the array. Because a node is subdivided before
// its children are visited, all descendants of a node are contiguous in the array, and children always come after
// their parent. The cells are the leaf nodes, numbered in the order in which they occur in the node array.
//
// The neighbor of a node across one of its walls is the node at the same level in the tree that touches the wall, if
// such a node exists. Otherwise, i.e. if the tree is less refined on the other side of the wall, it is the leaf node
// at a lower level that covers the other side of the wall. The neighbor is not defined for a wall on the boundary of
// the domain, in which case its index is -1. The neighbors of the two children of a node can be derived from the
// neighbors of the parent, so that the neighbors can be established with a single top-down pass through the array.
class BinTreeSpatialGrid::Tree
{
public:
    // this constructor creates the tree for the specified extent of the domain and the specified subdivision flags in
    // depth-first order, as returned by the readTopology() function
    Tree(const Box& extent, const vector<char>& flagv);

    // this function returns a pointer to the array of nodes; the first node is the root node
    const Node* nodes() const { return _nodev.data(); }

    // this function returns the number of cells in the tree, i.e. the number of leaf nodes
    int numCells() const { return static_cast<int>(_idv.size()); }

    // this function returns the leaf node corresponding to the cell with index m
    const Node& cellNode(int m) const { return _nodev[_idv[m]]; }

private:
    vector<Node> _nodev;  // all nodes in the tree; the first one is the root node
    vector<int> _idv;     // index in _nodev of the leaf node for each cell (i.e. for each cell index)
};

////////////////////////////////////////////////////////////////////

BinTreeSpatialGrid::Tree::Tree(const Box& extent, const vector<char>& flagv)
{
    // the node and cell indices are 32-bit integers
    if (flagv.size() >= static_cast<size_t>(std::numeric_limits<int>::max()))
        throw FATALERROR("Topology input file specifies too many nodes for a binary tree spatial grid");

    // create the nodes by subdividing depth-first, so that a subdivided node is followed by its two children
    _nodev.reserve(flagv.size());
    _nodev.emplace_back(extent, 0);
    vector<int> pendingv{0};  // indices of the nodes still waiting to be visited; the next one to visit is last
    for (char subdivide : flagv)
    {
        if (pendingv.empty()) throw FATALERROR("Topology input file has improper format and/or missing data");
        int n = pendingv.back();
        pendingv.pop_back();
        if (subdivide)
        {
            int c = static_cast<int>(_nodev.size());
            _nodev[n].setChild(c);
            _nodev.push_back(_nodev[n].lowerChild());
            _nodev.push_back(_nodev[n].upperChild());
            pendingv.push_back(c + 1);  // visit the lower child, and everything below it, before the upper child
            pendingv.push_back(c);
        }
    }
    if (!pendingv.empty()) throw FATALERROR("Topology input file has improper format and/or missing data");

    // establish the neighbors of the children of each nonleaf node, using the neighbors of the node itself (which
    // are known because the neighbors of the root node are all absent, and the nodes are in parent-first order)
    int numNodes = static_cast<int>(_nodev.size());
    for (int p = 0; p != numNodes; ++p)
    {
        const Node& parent = _nodev[p];
        if (parent.isLeaf()) continue;
        int c0 = parent.child();  // the lower child
        int c1 = c0 + 1;          // the upper child
        int a = parent.axis();

        // the children touch each other across the wall perpendicular to the split axis
        _nodev[c0].setNeighbor(2 * a + 1, c1);
        _nodev[c1].setNeighbor(2 * a, c0);

        // the other walls of the children coincide with a wall of the parent; both of the children touch the walls
        // perpendicular to the other axes, while the split axis has a single child at each of the parent's walls
        for (int w = 0; w != 6; ++w)
        {
            // if the parent has a neighbor at the same level that is subdivided, the neighbor at that wall is one of
            // the children of the parent's neighbor (which is split along the same axis, because it is at the same
            // level); otherwise it is the same node as for the parent
            int q = parent.neighbor(w);
            bool subdivided = q >= 0 && !_nodev[q].isLeaf() && _nodev[q].level() == parent.level();
            int q0 = subdivided ? _nodev[q].child() : q;      // lower child of the neighbor, or the neighbor itself
            int q1 = subdivided ? _nodev[q].child() + 1 : q;  // upper child of the neighbor, or the neighbor itself

            if (w / 2 != a)  // wall perpendicular to one of the other axes
            {
                _nodev[c0].setNeighbor(w, q0);
                _nodev[c1].setNeighbor(w, q1);
            }
            else if (w % 2 == 0)  // lower wall along the split axis: touched by the lower child only
            {
                _nodev[c0].setNeighbor(w, q1);  // adjacent to the upper child of the neighbor
            }
            else  // upper wall along the split axis: touched by the upper child only
            {
                _nodev[c1].setNeighbor(w, q0);  // adjacent to the lower child of the neighbor
            }
        }
    }

    // determine the cell indices and the leaf node for each cell
    for (int n = 0; n != numNodes; ++n)
    {
        if (_nodev[n].isLeaf())
        {
            _nodev[n].setCell(static_cast<int>(_idv.size()));
            _idv.push_back(n);
        }
    }
}

////////////////////////////////////////////////////////////////////

void BinTreeSpatialGrid::setupSelfAfter()
{
    BoxSpatialGrid::setupSelfAfter();

    // determine a small fraction relative to the spatial extent of the grid; used during path traversal
    _eps = 1e-12 * extent().diagonal();

    // open the input file
    string filepath = find<FilePaths>()->input(_filename);
    std::ifstream infile = System::ifstream(filepath);
    if (!infile) throw FATALERROR("Could not open the spatial tree grid topology text file " + filepath);

    // read the topology from the file
    Log* log = find<Log>();
    log->info("Constructing the binary tree spatial grid...");
    log->info("Reading tree topology from text file " + filepath + "...");
    vector<char> flagv = readTopology(infile);
    log->info("Done reading tree topology");

    // construct the tree
    _tree = std::make_shared<const Tree>(extent(), flagv);

    // determine the number of cells at each level in the tree hierarchy
    vector<int> countv;
    int numCells = _tree->numCells();
    for (int m = 0; m != numCells; ++m)
    {
        int level = _tree->cellNode(m).level();
        if (level + 1 > static_cast<int>(countv.size())) countv.resize(level + 1);
        countv[level]++;
    }

    // log these statistics, including a basic histogram
    log->info("Finished construction of the binary tree spatial grid");
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

int BinTreeSpatialGrid::numCells() const
{
    return _tree->numCells();
}

////////////////////////////////////////////////////////////////////

Box BinTreeSpatialGrid::cellBox(int m) const
{
    return _tree->cellNode(m).extent();
}

////////////////////////////////////////////////////////////////////

double BinTreeSpatialGrid::volume(int m) const
{
    return _tree->cellNode(m).extent().volume();
}

////////////////////////////////////////////////////////////////////

double BinTreeSpatialGrid::diagonal(int m) const
{
    return _tree->cellNode(m).extent().diagonal();
}

////////////////////////////////////////////////////////////////////

int BinTreeSpatialGrid::cellIndex(Position bfr) const
{
    const Node* nodes = _tree->nodes();
    int n = locate(nodes, bfr.x(), bfr.y(), bfr.z());
    return n >= 0 ? nodes[n].cell() : -1;
}

////////////////////////////////////////////////////////////////////

Position BinTreeSpatialGrid::centralPositionInCell(int m) const
{
    return Position(_tree->cellNode(m).extent().center());
}

////////////////////////////////////////////////////////////////////

Position BinTreeSpatialGrid::randomPositionInCell(int m) const
{
    return random()->position(_tree->cellNode(m).extent());
}

//////////////////////////////////////////////////////////////////////

class BinTreeSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const BinTreeSpatialGrid* _grid{nullptr};
    const Node* _nodes{nullptr};  // the array of nodes in the tree, cached for speed
    int _n{-1};                   // index of the leaf node containing the current position

    // the distance to the wall of a node that the path can cross along each axis is (wall - position) * inverse + bias;
    // these values depend only on the direction of the path, and are set up once for each path (see next())
    std::array<int, 3> _wallv{{1, 3, 5}};        // the wall (0-5) of a node through which the path can leave it
    std::array<double, 3> _invv{{0., 0., 0.}};   // the reciprocal of the direction component
    std::array<double, 3> _biasv{{0., 0., 0.}};  // zero, or a huge distance for a component that is zero

public:
    MySegmentGenerator(const BinTreeSpatialGrid* grid) : _grid(grid), _nodes(grid->_tree->nodes()) {}

    bool next() override
    {
        switch (state())
        {
            case State::Unknown:
            {
                // try moving the photon packet inside the grid; if this is impossible, return an empty path
                if (!moveInside(_grid->extent(), _grid->_eps)) return false;

                // get the node containing the current location;
                _n = locate(_nodes, rx(), ry(), rz());

                // The path can leave a node through the upper wall along each axis if the direction component is
                // positive, or through the lower wall otherwise, at a distance (wall - position) / component. To avoid
                // a division and any special cases in each step, we calculate the reciprocal of each component just
                // once for the path. A path parallel to an axis never crosses the walls perpendicular to that axis. We
                // represent this by a zero reciprocal and a huge offset, so that the distance to such a wall is a
                // constant that is never the shortest one. This avoids dividing by zero (which would yield a negative
                // infinity for a negative zero component and not-a-number for a position on the wall).
                const std::array<double, 3> kv{{kx(), ky(), kz()}};
                for (int a = 0; a != 3; ++a)
                {
                    bool nonzero = fabs(kv[a]) > 1e-300;
                    _wallv[a] = 2 * a + (kv[a] > 0. ? 1 : 0);
                    _invv[a] = nonzero ? 1. / kv[a] : 0.;
                    _biasv[a] = nonzero ? 0. : DBL_MAX;
                }

                // if the photon packet started outside the grid, return the corresponding nonzero-length segment;
                // otherwise fall through to determine the first actual segment
                if (ds() > 0.) return true;
            }

            // intentionally falls through
            case State::Inside:
            {
                // determine the segment from the current position to the first cell wall
                // and adjust the position and cell indices accordingly
                const Node& node = _nodes[_n];

                // The next node is one of the three nodes across the walls through which the path can leave this node.
                // Fetching a node from memory takes long compared to the calculations below, and the nodes of a large
                // tree are usually not in the cache. So we ask for all three candidates to be loaded right now, so that
                // this wait overlaps with the calculations. A wall on the boundary of the domain has no neighbor (-1),
                // in which case we harmlessly prefetch the root node instead.
                for (int a = 0; a != 3; ++a) prefetchNode(_nodes + std::max(node.neighbor(_wallv[a]), 0));

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
                propagater(ds + _grid->_eps);
                setSegment(node.cell(), ds);

                // find the new node by following the link across the crossed wall; this should not fail unless the new
                // location is outside the grid, however on rare occasions it fails due to rounding errors (e.g. in a
                // corner), in which case the function falls back to top-down search
                int oldn = _n;
                _n = locateAcross(_nodes, oldn, wall, rx(), ry(), rz());

                // if we're stuck in the same node,
                // try to escape by advancing the position to the next representable coordinates
                if (_n == oldn)
                {
                    propagateToNextAfter();
                    _n = locate(_nodes, rx(), ry(), rz());
                }

                // if we're outside the domain or still stuck in the same node, terminate the path
                if (_n < 0 || _n == oldn) setState(State::Outside);
                return true;
            }

            case State::Outside:
            {
            }
        }
        return false;
    }
};

////////////////////////////////////////////////////////////////////

std::unique_ptr<PathSegmentGenerator> BinTreeSpatialGrid::createPathSegmentGenerator() const
{
    return std::make_unique<MySegmentGenerator>(this);
}

////////////////////////////////////////////////////////////////////

void BinTreeSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
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

void BinTreeSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
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

void BinTreeSpatialGrid::write_yz(SpatialGridPlotFile* outfile) const
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

void BinTreeSpatialGrid::write_xyz(SpatialGridPlotFile* outfile) const
{
    // determine the number of cells at each level in the tree hierarchy
    vector<int> countv;
    int nCells = numCells();
    for (int m = 0; m != nCells; ++m)
    {
        int level = _tree->cellNode(m).level();
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
        const Node& node = _tree->cellNode(m);
        if (node.level() <= highestWriteLevel)
        {
            Box box = node.extent();
            outfile->writeCube(box.xmin(), box.ymin(), box.zmin(), box.xmax(), box.ymax(), box.zmax());
        }
    }
}

////////////////////////////////////////////////////////////////////

bool BinTreeSpatialGrid::offersInterface(const std::type_info& interfaceTypeInfo) const
{
    if (interfaceTypeInfo == typeid(DensityInCellInterface)) return BoxCellDensityMixIn::offersInterface();
    return BoxSpatialGrid::offersInterface(interfaceTypeInfo);
}

////////////////////////////////////////////////////////////////////
