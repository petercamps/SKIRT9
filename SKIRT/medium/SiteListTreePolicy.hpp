/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SITELISTTREEPOLICY_HPP
#define SITELISTTREEPOLICY_HPP

#include "Box.hpp"
#include "TreePolicy.hpp"
class TreeSpatialGrid;

//////////////////////////////////////////////////////////////////////

/** SiteListTreePolicy is a tree subdivision policy based on a list of site positions, such as the
    positions of the particles, sites or cells in an imported medium distribution. The site
    positions are loaded from a text column file if the \em filename property is nonempty;
    otherwise they are obtained from the first medium component in the medium system (in
    configuration order) that offers a site list, either itself (an imported medium) or through its
    geometry (a geometric medium with an imported geometry). See the ImportedMedium and
    ImportedGeometry classes and the various Snapshot subclasses for more information on the site
    positions returned by imported media and geometries.

    The input file, if any, contains a site position on each line, given as three columns with the
    \f$x\f$, \f$y\f$, and \f$z\f$ coordinates. The default unit is pc; this can be overridden by
    providing column information in the header of the file, as described for the TextInFile class.
    Sites outside of the spatial domain of the grid are ignored.

    The policy requests subdivision in two steps. In a first step, the tree is subdivided in such a
    way that each leaf node contains at most one of the sites in the list. Subsequently, each of
    these leaf nodes, including the empty ones, is further subdivided a fixed number of times, as
    configured by the \em numExtraLevels property. The minimum and maximum subdivision levels of
    the grid override these criteria. Nodes are always subdivided up to the minimum level, so that
    the first step never produces leaf nodes coarser than that level, and the extra levels are
    counted from the minimum level for a site that is isolated at a lower level. Nodes are never
    subdivided beyond the maximum level, so that a leaf node at that level may contain multiple
    sites. When this policy is the only policy of the grid, the resulting tree is purely geometric,
    i.e. it does not depend on the random number sequence.

    To evaluate a node independently of any other node, the policy constructs a reference tree
    during setup. This tree mirrors the first step described above, i.e. it is subdivided only
    where a node contains more than one site, using the same subdivision rule as the grid. While
    constructing the grid, the policy descends the reference tree from its root toward the center
    of the node being evaluated, up to the level of that node. If the reference tree is still
    subdivided at that level, the node contains multiple sites and must be subdivided. Otherwise,
    the descent ended at a leaf of the reference tree, and the node is subdivided if its level is
    below the level of that leaf plus the number of extra levels. Because the reference tree is not
    subdivided where a node holds at most one site, its leaves may be coarser than the minimum
    level of the grid, while the grid always subdivides down to that minimum level. In that case,
    the extra levels are counted from the minimum level instead. */
class SiteListTreePolicy : public TreePolicy
{
    ITEM_CONCRETE(SiteListTreePolicy, TreePolicy, "a tree subdivision policy using a list of site positions")
        ATTRIBUTE_TYPE_DISPLAYED_IF(SiteListTreePolicy, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the site positions (empty to use an imported medium)")
        ATTRIBUTE_REQUIRED_IF(filename, "!SiteListInterface")

        PROPERTY_INT(numExtraLevels, "the number of additional subdivision levels")
        ATTRIBUTE_MIN_VALUE(numExtraLevels, "0")
        ATTRIBUTE_MAX_VALUE(numExtraLevels, "30")
        ATTRIBUTE_DEFAULT_VALUE(numExtraLevels, "0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function obtains the site positions and constructs the reference tree as described in
        the class header. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns true if the specified node needs to be subdivided according to the
        criteria described in the class header, and false otherwise. */
    bool needsSubdivide(TreeNodeEvaluation& node) const override;

    //======================== Data Members ========================

private:
    // data members initialized by setupSelfBefore()
    const TreeSpatialGrid* _grid{nullptr};  // the tree grid, which offers its subdivision rule
    Box _extent;                            // the spatial domain of the grid
    int _minLevel{0};                       // the minimum level of the grid
    vector<int> _childv;  // index of the first child for each node in the reference tree, or -1 for a leaf
};

//////////////////////////////////////////////////////////////////////

#endif
