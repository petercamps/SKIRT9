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

    Locating the cell containing a given position and generating path segments are implemented by
    the TreeSpatialGrid base class for both tree types. */
class OctTreeSpatialGrid : public TreeSpatialGrid
{
    ITEM_CONCRETE(OctTreeSpatialGrid, TreeSpatialGrid, "an octtree spatial grid (8 children per node)")
    ITEM_END()

    //=========== Functions implementing the tree type ===========

public:
    /** This function returns the number of children of a nonleaf node, i.e. 8. */
    int numChildren() const override;

    /** This function returns the extent of the child with index \f$c\f$ of a node with the
        specified extent and level, as described in the class header. */
    Box childExtent(const Box& extent, int level, int c) const override;

    /** This function returns the index of the child that contains the specified position, for a
        node with the specified extent and level. A position on a splitting plane is assigned to the
        upper child. */
    int childIndex(const Box& extent, int level, Vec position) const override;

protected:
    /** This function establishes the neighbor links for all nodes in the specified array, as
        described in the class header. */
    void linkNeighbors(vector<Node>& nodes) const override;
};

//////////////////////////////////////////////////////////////////////

#endif
