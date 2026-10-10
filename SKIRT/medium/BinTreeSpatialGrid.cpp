/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BinTreeSpatialGrid.hpp"

////////////////////////////////////////////////////////////////////

int BinTreeSpatialGrid::numChildren() const
{
    return 2;
}

////////////////////////////////////////////////////////////////////

Box BinTreeSpatialGrid::childExtent(const Box& extent, int level, int c) const
{
    // the lower child (c=0) extends from the lower wall of the parent to the splitting plane perpendicular to the
    // axis for this level, and the upper child (c=1) from the splitting plane to the upper wall of the parent
    Node node(extent, level);
    int a = node.axis();
    node.setWall(c ? 2 * a : 2 * a + 1, node.center(a));
    return node.extent();
}

////////////////////////////////////////////////////////////////////

int BinTreeSpatialGrid::childIndex(const Box& extent, int level, Vec position) const
{
    Node node(extent, level);
    int a = node.axis();
    const double r[3] = {position.x(), position.y(), position.z()};
    return r[a] < node.center(a) ? 0 : 1;
}

////////////////////////////////////////////////////////////////////

void BinTreeSpatialGrid::linkNeighbors(vector<Node>& nodes) const
{
    // establish the neighbors of the children of each nonleaf node, using the neighbors of the node itself (which
    // are known because the neighbors of the root node are all absent, and the nodes are in parent-first order)
    int numNodes = static_cast<int>(nodes.size());
    for (int p = 0; p != numNodes; ++p)
    {
        const Node& parent = nodes[p];
        if (parent.isLeaf()) continue;
        int c0 = parent.child();  // the lower child
        int c1 = c0 + 1;          // the upper child
        int a = parent.axis();

        // the children touch each other across the wall perpendicular to the split axis
        nodes[c0].setNeighbor(2 * a + 1, c1);
        nodes[c1].setNeighbor(2 * a, c0);

        // the other walls of the children coincide with a wall of the parent; both of the children touch the walls
        // perpendicular to the other axes, while the split axis has a single child at each of the parent's walls
        for (int w = 0; w != 6; ++w)
        {
            // if the parent has a neighbor at the same level that is subdivided, the neighbor at that wall is one of
            // the children of the parent's neighbor (which is split along the same axis, because it is at the same
            // level); otherwise it is the same node as for the parent
            int q = parent.neighbor(w);
            bool subdivided = q >= 0 && !nodes[q].isLeaf() && nodes[q].level() == parent.level();
            int q0 = subdivided ? nodes[q].child() : q;      // lower child of the neighbor, or the neighbor itself
            int q1 = subdivided ? nodes[q].child() + 1 : q;  // upper child of the neighbor, or the neighbor itself

            if (w / 2 != a)  // wall perpendicular to one of the other axes
            {
                nodes[c0].setNeighbor(w, q0);
                nodes[c1].setNeighbor(w, q1);
            }
            else if (w % 2 == 0)  // lower wall along the split axis: touched by the lower child only
            {
                nodes[c0].setNeighbor(w, q1);  // adjacent to the upper child of the neighbor
            }
            else  // upper wall along the split axis: touched by the upper child only
            {
                nodes[c1].setNeighbor(w, q0);  // adjacent to the lower child of the neighbor
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
