/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TREEPOLICY_HPP
#define TREEPOLICY_HPP

#include "SimulationItem.hpp"
class TreeNodeEvaluation;

//////////////////////////////////////////////////////////////////////

/** TreePolicy is an abstract class that represents a criterion for subdividing the nodes of a
    spatial tree grid, as a service to the TreeSpatialGrid class. A tree grid holds a list of
    policies. While constructing the tree, the grid asks each policy in turn whether a given node
    needs to be subdivided, and subdivides the node as soon as one of the policies asks for it. As
    a result, the policies can be combined freely. The minimum and maximum subdivision levels
    configured for the grid override the outcome of the policies.

    The public interface of a policy consists of a single function, needsSubdivide(), that
    evaluates the criterion for a given node. The node is represented by a TreeNodeEvaluation
    instance, which offers its spatial extent and level in the tree, and the properties of the
    media in the node. These properties are calculated when first requested and then cached, so
    that the policies evaluating the same node share them, including the random positions at
    which the densities are sampled. The criterion for a node depends only on the node itself, not
    on any other nodes. The tree grid evaluates the nodes at a given level in parallel, each thread
    with its own TreeNodeEvaluation instance, so the needsSubdivide() function must be
    thread-safe. */
class TreePolicy : public SimulationItem
{
    ITEM_ABSTRACT(TreePolicy, SimulationItem, "a spatial tree grid subdivision policy")
    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function must be implemented in a subclass. It returns true if the specified node
        needs to be subdivided according to this policy, and false otherwise. The function may be
        called in parallel from multiple threads (each with its own TreeNodeEvaluation instance),
        and must thus be thread-safe. */
    virtual bool needsSubdivide(TreeNodeEvaluation& node) const = 0;
};

//////////////////////////////////////////////////////////////////////

#endif
