/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CELLTREEPOLICY_HPP
#define CELLTREEPOLICY_HPP

#include "TreePolicy.hpp"
class CellSnapshot;

//////////////////////////////////////////////////////////////////////

/** CellTreePolicy is a tree subdivision policy that ensures that the grid is nowhere coarser
    than an imported list of cuboidal cells lined up with the coordinate axes. A node is subdivided
    if it is wider, in any of the spatial directions, than the cell that contains the center of
    the node. If multiple cells contain the center, the cell listed first in the imported file is
    used, as for the other properties of a cell snapshot (see the CellSnapshot class). Nodes whose
    center is not inside any of the cells are not subdivided by this policy.
    The criterion is approximate, because only the cell at the center of each node is considered;
    a node may thus remain coarser than smaller cells elsewhere within its extent.

    If the \em filename property is nonempty, the policy imports its own cells from the specified
    file, in the format described for the CellSnapshot class: each line specifies the coordinates
    \f$x_\text{min}\f$, \f$y_\text{min}\f$, \f$z_\text{min}\f$, \f$x_\text{max}\f$,
    \f$y_\text{max}\f$, and \f$z_\text{max}\f$ of a cell (any further columns are ignored).
    Otherwise, the policy uses the cells imported by the first medium component in the medium
    system (in configuration order) that offers them, either itself (a cell medium) or through its
    geometry (a geometric medium with a cell geometry). This avoids a second copy of the same
    cells, and guarantees that the policy uses the cells of the actual medium. */
class CellTreePolicy : public TreePolicy
{
    ITEM_CONCRETE(CellTreePolicy, TreePolicy, "a tree subdivision policy following a list of cuboidal cells")
        ATTRIBUTE_TYPE_DISPLAYED_IF(CellTreePolicy, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the cells (empty to use an imported medium)")
        ATTRIBUTE_REQUIRED_IF(filename, "!CellMeshInterface")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor deletes the imported cells, if any. */
    ~CellTreePolicy();

protected:
    /** This function imports the cells or locates the cells imported by a medium component. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns true if the specified node is wider than the cell containing its
        center in any of the spatial directions, and false otherwise. */
    bool needsSubdivide(TreeNodeEvaluation& node) const override;

    //======================== Data Members ========================

private:
    // data members initialized by setupSelfBefore()
    CellSnapshot* _ownMesh{nullptr};     // the cells imported by this policy, if any
    const CellSnapshot* _mesh{nullptr};  // the cells used by this policy
};

//////////////////////////////////////////////////////////////////////

#endif
