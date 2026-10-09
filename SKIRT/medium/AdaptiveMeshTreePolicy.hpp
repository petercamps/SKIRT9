/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADAPTIVEMESHTREEPOLICY_HPP
#define ADAPTIVEMESHTREEPOLICY_HPP

#include "TreePolicy.hpp"
class AdaptiveMeshSnapshot;

//////////////////////////////////////////////////////////////////////

/** AdaptiveMeshTreePolicy is a tree subdivision policy that ensures that the grid is nowhere
    coarser than an imported adaptive mesh. A node is subdivided if it is wider, in any of the
    spatial directions, than the leaf cell of the adaptive mesh that contains the center of the
    node. For an adaptive mesh in which each nonleaf node is split in two along each axis, with a
    domain that coincides with the domain of an octtree grid, the nodes of the grid line up with
    the cells of the mesh, so that the grid exactly reproduces the mesh where it is not refined
    further by other policies. For other meshes, the criterion is approximate, because only the
    mesh cell at the center of each node is considered.

    If the \em filename property is nonempty, the policy imports its own adaptive mesh from the
    specified file, in the format described for the AdaptiveMeshSnapshot class (data values for
    the leaf cells are ignored), using the spatial domain of the grid as the domain of the mesh.
    Otherwise, the policy uses the adaptive mesh imported by the first medium component in the
    medium system (in configuration order) that offers one, either itself (an adaptive mesh
    medium) or through its geometry (a geometric medium with an adaptive mesh geometry). This
    avoids a second copy of the same mesh, and guarantees that the policy uses the mesh of the
    actual medium. */
class AdaptiveMeshTreePolicy : public TreePolicy
{
    ITEM_CONCRETE(AdaptiveMeshTreePolicy, TreePolicy, "a tree subdivision policy following an adaptive mesh (AMR grid)")
        ATTRIBUTE_TYPE_DISPLAYED_IF(AdaptiveMeshTreePolicy, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the adaptive mesh (empty to use an imported medium)")
        ATTRIBUTE_REQUIRED_IF(filename, "!AdaptiveMeshInterface")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor deletes the imported adaptive mesh, if any. */
    ~AdaptiveMeshTreePolicy();

protected:
    /** This function imports the adaptive mesh, or locates the adaptive mesh imported by a medium
        component. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns true if the specified node is wider than the cell containing its
        center in any of the spatial directions, and false otherwise. */
    bool needsSubdivide(TreeNodeEvaluation& node) const override;

    //======================== Data Members ========================

private:
    // data members initialized by setupSelfBefore()
    AdaptiveMeshSnapshot* _ownMesh{nullptr};     // the adaptive mesh imported by this policy, if any
    const AdaptiveMeshSnapshot* _mesh{nullptr};  // the adaptive mesh used by this policy
};

//////////////////////////////////////////////////////////////////////

#endif
