/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OctTreeSpatialGrid.hpp"
#include "PathSegmentGenerator.hpp"
#include "Prefetch.hpp"

////////////////////////////////////////////////////////////////////

int OctTreeSpatialGrid::numChildren() const
{
    return 8;
}

////////////////////////////////////////////////////////////////////

Box OctTreeSpatialGrid::childExtent(const Box& extent, int level, int c) const
{
    // the child with index c covers the lower or upper half of the parent along axis a depending on bit a of c
    const Node parent(extent, level);
    Node child(extent, level + 1);
    for (int a = 0; a != 3; ++a)
    {
        if ((c & (1 << a)) != 0)
            child.setWall(2 * a, parent.center(a));  // upper half: the lower wall moves to the center
        else
            child.setWall(2 * a + 1, parent.center(a));  // lower half: the upper wall moves to the center
    }
    return child.extent();
}

////////////////////////////////////////////////////////////////////

int OctTreeSpatialGrid::childIndex(const Box& extent, int /*level*/, Vec position) const
{
    Vec center = extent.center();
    return (position.x() < center.x() ? 0 : 1) + (position.y() < center.y() ? 0 : 2)
           + (position.z() < center.z() ? 0 : 4);
}

////////////////////////////////////////////////////////////////////

void OctTreeSpatialGrid::linkNeighbors(vector<Node>& nodes) const
{
    // establish the neighbors of the children of each nonleaf node, using the neighbors of the node itself (which
    // are known because the neighbors of the root node are all absent, and the nodes are in parent-first order)
    int numNodes = static_cast<int>(nodes.size());
    for (int p = 0; p != numNodes; ++p)
    {
        const Node& parent = nodes[p];
        if (parent.isLeaf()) continue;
        int c0 = parent.child();  // the first child

        for (int i = 0; i != 8; ++i)
        {
            Node& child = nodes[c0 + i];
            for (int a = 0; a != 3; ++a)
            {
                int bit = 1 << a;
                bool upper = (i & bit) != 0;  // true if the child covers the upper half of the parent along axis a

                // the wall on the side of the parent's center touches the sibling on the other side of the center
                child.setNeighbor(upper ? 2 * a : 2 * a + 1, c0 + (i ^ bit));

                // the wall on the other side coincides with a wall of the parent; if the parent has a neighbor at the
                // same level that is subdivided, the neighbor at that wall is the child of the parent's neighbor in
                // the mirrored octant across that wall; otherwise it is the same node as for the parent
                int w = upper ? 2 * a + 1 : 2 * a;
                int q = parent.neighbor(w);
                bool subdivided = q >= 0 && !nodes[q].isLeaf() && nodes[q].level() == parent.level();
                child.setNeighbor(w, subdivided ? nodes[q].child() + (i ^ bit) : q);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    using Node = TreeSpatialGrid::Node;

    // This function returns the index of the leaf node that contains the specified position, descending the tree from
    // the node with index n in the specified array of nodes, which should contain the position. A position on a
    // splitting plane is assigned to the child on the upper side of the plane. If the position is not inside the
    // starting node, the result is not necessarily the leaf node closest to the position, and the caller must verify
    // the result.
    int descend(const vector<Node>& nodes, int n, double x, double y, double z)
    {
        while (!nodes[n].isLeaf())
        {
            const Node& node = nodes[n];
            n = node.child() + (x < node.center(0) ? 0 : 1) + (y < node.center(1) ? 0 : 2)
                + (z < node.center(2) ? 0 : 4);
        }
        return n;
    }

    // This function returns the index of the leaf node that contains the specified position by searching the tree
    // top-down from its root node, or -1 if the position is outside of the domain.
    int locate(const vector<Node>& nodes, double x, double y, double z)
    {
        if (!nodes[0].contains(x, y, z)) return -1;
        return descend(nodes, 0, x, y, z);
    }

    // This function returns the index of the leaf node that contains the specified position, which must be located
    // just across the specified wall (0-5) of the leaf node with index n. It proceeds from the neighbor across that
    // wall, and returns -1 if the position is outside of the domain. If the neighbor does not contain the position due
    // to numerical inaccuracies (for example, when a path passes very close to the edge or corner of a node), it uses
    // a top-down search as a fall-back.
    int locateAcross(const vector<Node>& nodes, int n, int wall, double x, double y, double z)
    {
        int m = nodes[n].neighbor(wall);
        if (m >= 0)
        {
            m = descend(nodes, m, x, y, z);
            if (nodes[m].contains(x, y, z)) return m;
        }
        return locate(nodes, x, y, z);
    }

    // This class implements the path segment generator for an octtree grid, as described in the class header.
    class SegmentGenerator : public PathSegmentGenerator
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
                    // if the top-down search would end in the leaf node for the known initial cell, start from that
                    // node without searching, and determine the first segment
                    int n = _cellNodes[initialCellIndex()];
                    if (_nodes[n].owns(rx(), ry(), rz()))
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
                    _n = locate(_nodes, rx(), ry(), rz());
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
    };
}

////////////////////////////////////////////////////////////////////

int OctTreeSpatialGrid::cellIndex(Position bfr) const
{
    int n = locate(nodes(), bfr.x(), bfr.y(), bfr.z());
    return n >= 0 ? nodes()[n].cell() : -1;
}

////////////////////////////////////////////////////////////////////

std::unique_ptr<PathSegmentGenerator> OctTreeSpatialGrid::createPathSegmentGenerator() const
{
    return std::make_unique<SegmentGenerator>(extent(), eps(), nodes(), cellNodeIndices());
}

////////////////////////////////////////////////////////////////////
