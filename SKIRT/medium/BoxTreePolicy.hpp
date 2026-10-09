/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BOXTREEPOLICY_HPP
#define BOXTREEPOLICY_HPP

#include "Box.hpp"
#include "TreePolicy.hpp"

//////////////////////////////////////////////////////////////////////

/** BoxTreePolicy is a tree subdivision policy that applies a nested policy only to the nodes that
    intersect a given box. This can be used, for example, to specify a higher resolution in a given
    region of interest. The nested policy can be any tree subdivision policy, including another
    BoxTreePolicy instance.

    A node is subdivided if it intersects the box (borders included) and the nested policy asks for
    subdivision. A node intersecting the box is thus subdivided according to the criteria of the
    nested policy, even if it is mostly outside of the box. For nodes that do not intersect the
    box, this policy never asks for subdivision. To obtain a coarser resolution outside the box,
    the grid lists one or more other policies that cover the full domain alongside this one, with
    less strict criteria than the nested policy. Because a node is subdivided as soon as one of the
    policies asks for it, the stricter criteria of the nested policy then dominate inside the box.
    Several BoxTreePolicy instances can be listed to refine multiple regions, each with its own
    criteria. The minimum and maximum levels of the grid apply to the whole domain.

    It is not meaningful to specify a box that extends outside of the spatial grid domain, because
    none of the tree nodes will intersect those outside areas. Therefore, a warning is issued if the
    box is not fully inside the spatial grid domain. */
class BoxTreePolicy : public TreePolicy
{
    ITEM_CONCRETE(BoxTreePolicy, TreePolicy, "a tree subdivision policy applying a nested policy inside a box")
        ATTRIBUTE_TYPE_DISPLAYED_IF(BoxTreePolicy, "Level2")

        PROPERTY_DOUBLE(minX, "the start point of the box in the X direction")
        ATTRIBUTE_QUANTITY(minX, "length")

        PROPERTY_DOUBLE(maxX, "the end point of the box in the X direction")
        ATTRIBUTE_QUANTITY(maxX, "length")

        PROPERTY_DOUBLE(minY, "the start point of the box in the Y direction")
        ATTRIBUTE_QUANTITY(minY, "length")

        PROPERTY_DOUBLE(maxY, "the end point of the box in the Y direction")
        ATTRIBUTE_QUANTITY(maxY, "length")

        PROPERTY_DOUBLE(minZ, "the start point of the box in the Z direction")
        ATTRIBUTE_QUANTITY(minZ, "length")

        PROPERTY_DOUBLE(maxZ, "the end point of the box in the Z direction")
        ATTRIBUTE_QUANTITY(maxZ, "length")

        PROPERTY_ITEM(policy, TreePolicy, "the policy applied to nodes intersecting the box")
        ATTRIBUTE_DEFAULT_VALUE(policy, "DensityTreePolicy")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the box is not empty, and warns if it is not inside the spatial
        domain of the grid. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns true if the specified node intersects the box and the nested policy
        asks for subdivision, and false otherwise. */
    bool needsSubdivide(TreeNodeEvaluation& node) const override;

    //======================== Data Members ========================

private:
    // data members initialized by setupSelfBefore()
    Box _box;
};

//////////////////////////////////////////////////////////////////////

#endif
