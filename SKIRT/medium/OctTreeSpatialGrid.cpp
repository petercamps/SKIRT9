/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "OctTreeSpatialGrid.hpp"

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
