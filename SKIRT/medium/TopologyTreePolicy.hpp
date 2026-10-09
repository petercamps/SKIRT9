/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TOPOLOGYTREEPOLICY_HPP
#define TOPOLOGYTREEPOLICY_HPP

#include "Box.hpp"
#include "TreePolicy.hpp"
class TreeSpatialGrid;

//////////////////////////////////////////////////////////////////////

/** TopologyTreePolicy is a tree subdivision policy that reproduces a tree topology loaded from a
    data file. In the most common use case, the data file has been created in a previous simulation
    through the TreeSpatialGridTopologyProbe. When this policy is the only policy of the grid, the
    tree grid in the current simulation will then be identical to the one in the previous
    simulation, provided that the grid is of the same type (binary tree or octtree) and that its
    range of subdivision levels is sufficiently wide: the minimum level should be zero (or at most
    the level of the coarsest leaf in the recorded tree), and the maximum level at least the depth
    of the recorded tree. The policy logs a warning if the maximum level is smaller than that
    depth.

    Reproducing a recorded tree can be useful in situations where multiple simulations are being
    performed on input models with the same (or a very similar) spatial distribution of the
    transfer medium. Constructing a tree grid based on the medium density distribution takes time
    (because of the large number of density samples required), and the tree is likely to differ
    slightly between various runs (because the density is sampled at random positions, and the
    random number sequence varies because of parallelization). Loading the tree from the topology
    data file is much faster and guarantees that all simulations use identical tree grids. A
    recorded tree can also be refined by configuring other policies alongside this one, because a
    node is subdivided as soon as one of the policies asks for it.

    After a brief descriptive header (lines starting with a # character), the input file contains
    lines with just a single integer number. The first line specifies the number of children for
    each nonleaf node (2 for a binary tree, 8 for an octtree, or 0 if the root node is not
    subdivided). The second line contains 1 if the root node is subdivided, or 0 if not. The
    following lines similarly contain 1 or 0 indicating subdivision for any children of the
    preceding node, recursively, in a depth-first traversal of the tree. The topology data is
    scale-free; the spatial domain is defined by the grid in the current simulation. If the number
    of children specified in the file differs from that of the grid (and is not zero), setup reports
    a fatal error.

    During setup, the policy loads the recorded topology into a reference tree. To evaluate a node
    of the grid being constructed, it descends the reference tree from its root toward the center of
    the node, using the subdivision rule of the grid, up to the level of the node. The node is
    subdivided if the reference tree reaches that level and is subdivided there; it is not
    subdivided if the reference node at that level is a leaf, or if the reference tree is less deep
    at that location. */
class TopologyTreePolicy : public TreePolicy
{
    ITEM_CONCRETE(TopologyTreePolicy, TreePolicy, "a tree subdivision policy using a topology loaded from file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TopologyTreePolicy, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the tree topology data")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function loads the recorded topology into a reference tree as described in the class
        header. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns true if the reference tree is subdivided at the location and level of
        the specified node, and false otherwise. */
    bool needsSubdivide(TreeNodeEvaluation& node) const override;

    //======================== Data Members ========================

private:
    // data members initialized by setupSelfBefore()
    const TreeSpatialGrid* _grid{nullptr};  // the tree grid, which offers its subdivision rule
    Box _extent;                            // the spatial domain of the grid
    vector<int> _childv;  // index of the first child for each node in the reference tree, or -1 for a leaf
};

//////////////////////////////////////////////////////////////////////

#endif
