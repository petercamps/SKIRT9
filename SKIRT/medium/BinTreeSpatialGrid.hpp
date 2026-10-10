/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BINTREESPATIALGRID_HPP
#define BINTREESPATIALGRID_HPP

#include "TreeSpatialGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** BinTreeSpatialGrid is a concrete subclass of the TreeSpatialGrid class representing a spatial
    grid with cuboidal cells organized in a binary tree. Each nonleaf node is divided into two equal
    halves along a plane perpendicular to one of the coordinate axes, alternating between the x, y
    and z axes when descending the tree: a node at level \f$l\f$ is split perpendicular to axis
    \f$l \bmod 3\f$. The lower child (the one with the lowest coordinates along the split axis)
    comes first. The tree construction, according to the configured subdivision policies, is
    described in the TreeSpatialGrid class.

    A binary tree grid can follow the spatial structure of the medium more closely than an octtree
    grid, because a node is split in just two halves. As a result, it usually needs fewer cells for
    the same resolution, but the tree is deeper.

    <b>Neighbors</b>

    The neighbor across a wall of a node is the node at the same level in the tree that touches the
    wall if there is one; otherwise it is the leaf node (at a lower level) that covers the other
    side of the wall. The neighbors of the two children of a node can be derived from the neighbors
    of the node itself, so that the neighbors are established in a single top-down pass through the
    tree.

    <b>Path segment generation</b>

    Locating the cell containing a given position and generating path segments are implemented by
    the TreeSpatialGrid base class for both tree types. */
class BinTreeSpatialGrid : public TreeSpatialGrid
{
    ITEM_CONCRETE(BinTreeSpatialGrid, TreeSpatialGrid, "a binary tree spatial grid (2 children per node)")
        ATTRIBUTE_TYPE_DISPLAYED_IF(BinTreeSpatialGrid, "Level2")
    ITEM_END()

    //=========== Functions implementing the tree type ===========

public:
    /** This function returns the number of children of a nonleaf node, i.e. 2. */
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
