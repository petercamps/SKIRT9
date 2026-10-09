/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef OCTTREESPATIALGRID_HPP
#define OCTTREESPATIALGRID_HPP

#include "TreeSpatialGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** OctTreeSpatialGrid is a concrete subclass of the TreeSpatialGrid class representing a spatial
    grid with cuboidal cells organized in an octtree. Each nonleaf node is divided into eight equal
    octants by three planes through its center, perpendicular to the coordinate axes. The child with
    index \f$i = i_x + 2 i_y + 4 i_z\f$ relative to the first child covers the lower half of the
    node along axis \f$a\f$ if \f$i_a=0\f$ and the upper half if \f$i_a=1\f$. The tree construction,
    according to the configured subdivision policies, is described in the TreeSpatialGrid class.

    <b>Neighbors</b>

    The neighbor across a wall of a node is the node at the same level in the tree that touches the
    wall if there is one; otherwise it is the leaf node (at a lower level) that covers the other
    side of the wall. The neighbors of the eight children of a node can be derived from the
    neighbors of the node itself, so that the neighbors are established in a single top-down pass
    through the tree.

    <b>Path segment generation</b>

    The path segment generator determines the cell that contains the starting position, and
    calculates the first wall of the cell that will be crossed. The path length \f$\Delta s\f$ is
    determined and the current position is moved to a new position along this path, a tiny fraction
    further than \f$\Delta s\f$, so that the new position is within the next cell. To determine that
    next cell, the generator follows the link to the neighbor across the crossed wall. If the
    neighbor turns out to have children of its own (because the tree is more refined on the other
    side of the wall), the generator descends into the neighbor until it reaches the leaf containing
    the new position. A top-down search starting at the root node is used to determine the initial
    cell and as a fall-back in rare cases where numerical inaccuracies would otherwise result in an
    inconsistent state.

    The quantities that depend only on the direction of the path, i.e. the reciprocal of each
    direction component and the walls that the path can cross, are calculated just once for each
    path. The nodes of a large tree are usually not in the processor cache, so that waiting for the
    next node to arrive from memory dominates the cost of each step. To hide this latency, the
    generator asks the processor to prefetch the (up to) three nodes that the path can move to next
    while it is still working on the current node (see the Prefetch namespace).

    A path that runs exactly along a cell boundary (for example, a path parallel to a coordinate
    axis through a position on a splitting plane) is ambiguous: the cells on either side of the
    boundary are equally valid choices. In such cases, the generator consistently selects the cell
    on the upper side of a splitting plane, as it does when locating the cell that contains a given
    position. */
class OctTreeSpatialGrid : public TreeSpatialGrid
{
    ITEM_CONCRETE(OctTreeSpatialGrid, TreeSpatialGrid, "an octtree spatial grid (8 children per node)")
    ITEM_END()

    //=========== Functions implementing the tree type ===========

public:
    /** This function returns the number of children of a nonleaf node, i.e. 8. */
    int numChildren() const override;

protected:
    /** This function appends the eight children of the specified node to the end of the specified
        list of nodes, as described in the class header. */
    void appendChildren(const Node& parent, std::deque<Node>& nodes) const override;

    /** This function establishes the neighbor links for all nodes in the specified array, as
        described in the class header. */
    void linkNeighbors(vector<Node>& nodes) const override;

    //======================== Other Functions =======================

public:
    /** This function returns the index of the cell that contains the position \f${\bf{r}}\f$. The
        search algorithm starts at the root node and selects the child node that contains the
        position. This procedure is repeated until the node is childless, i.e. until it is a leaf
        node that corresponds to an actual spatial cell. A position on a splitting plane is
        assigned to the child on the upper side of the plane. */
    int cellIndex(Position bfr) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for this grid, implemented as a
        PathSegmentGenerator subclass local to the implementation file. The algorithm is described
        in the class header. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;
};

//////////////////////////////////////////////////////////////////////

#endif
