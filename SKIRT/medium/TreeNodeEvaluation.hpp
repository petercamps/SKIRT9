/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TREENODEEVALUATION_HPP
#define TREENODEEVALUATION_HPP

#include "Box.hpp"
#include "MaterialMix.hpp"
#include "Position.hpp"
#include <array>
class MassInBoxInterface;
class Medium;
class MediumSystem;
class Random;

//////////////////////////////////////////////////////////////////////

/** A TreeNodeEvaluation instance represents a node of a spatial tree grid while the subdivision
    policies of the grid evaluate whether the node needs to be subdivided (see the TreePolicy and
    TreeSpatialGrid classes). It offers the extent and level of the node, and the properties of the
    media in the node that are needed by the policies.

    The properties of the media are calculated only when a policy requests them, and they are
    cached, so that the policies evaluating the same node share them. This includes the random
    positions at which the densities are sampled: they are drawn when the first density sample is
    requested, and all medium components and policies use the same sample positions. As a result,
    the evaluation of a node costs no more than necessary, regardless of the number of policies.

    The quantities for a dust medium component are expressed in terms of mass, and those for an
    electron or gas medium component in terms of number of entities. The amount of material of a
    medium component in the node is calculated exactly if the component offers the
    MassInBoxInterface. Otherwise, it is estimated from the density sampled in
    \f$N\f$ random positions \f$\bf{r}_i\f$ distributed uniformly across the volume \f$V\f$ of the
    node, \f[ M = \frac{V}{N} \sum_{i=1}^{N}\rho(\bf{r}_i). \f] The number of samples \f$N\f$ is
    configured in the SamplingOptions of the medium system.

    A tree grid evaluates the nodes at a given level in parallel. Each thread uses its own
    TreeNodeEvaluation instance, which it resets for each node, so that the cached information does
    not need to be protected against concurrent access, and the memory used for caching can be
    reused for the next node. */
class TreeNodeEvaluation
{
public:
    /** This constructor prepares an instance for evaluating nodes in the context of the specified
        medium system (which may be the null pointer if there is none), using the specified random
        number generator and number of density samples per node. The medium system must have been
        set up at least partially, so that its lists of medium component indices per material type
        are available. */
    TreeNodeEvaluation(const MediumSystem* ms, Random* random, int numSamples);

    /** This function prepares the instance for evaluating the node with the specified extent and
        level in the tree, discarding any information cached for the previous node. */
    void reset(const Box& box, int level);

    /** This function returns the extent of the node. */
    const Box& box() const { return _box; }

    /** This function returns the level of the node in the tree, with level zero for the root
        node. */
    int level() const { return _level; }

    /** This function returns the average density of the medium component with index \f$h\f$ in
        the node, i.e. the amount of material of the component in the node divided by its volume:
        the mass density for a dust component, or the number density for an electron or gas
        component. */
    double density(int h);

    /** This function returns the amount of material of the specified type in the node, summed over
        the medium components of that type: the mass for dust, or the number of entities for
        electrons and gas. For dust, mass is the appropriate quantity in case multiple dust medium
        components have a different mass per hydrogen atom. */
    double amount(MaterialMix::MaterialType type);

    /** This function returns the density of the material of the specified type at each of the
        sample positions in the node, summed over the medium components of that type: the mass
        density for dust, or the number density for electrons and gas. */
    const vector<double>& densitySamples(MaterialMix::MaterialType type);

    /** This function returns the density of the medium component with index \f$h\f$ at each of
        the sample positions in the node: the mass density for a dust component, or the number
        density for an electron or gas component. */
    const vector<double>& componentDensitySamples(int h);

    /** This function returns the random positions, distributed uniformly across the volume of the
        node, at which the densities are sampled. The positions are drawn when first requested,
        either through this function or through one of the functions returning density samples.
        A policy can use them to sample other quantities at the same positions. */
    const vector<Position>& samplePositions();

    /** This function returns a lowercase name for the specified material type ("dust",
        "electrons", or "gas"), for use in messages. */
    static string materialTypeName(MaterialMix::MaterialType type);

private:
    // returns the list of medium component indices for the specified material type
    const vector<int>& indices(MaterialMix::MaterialType type) const;

    // returns the amount of material of the specified medium component in the node
    double componentAmount(int h);

    // data members initialized by the constructor
    const MediumSystem* _ms{nullptr};
    Random* _random{nullptr};
    int _numSamples{0};
    vector<Medium*> _media;             // the medium components
    vector<MassInBoxInterface*> _mibv;  // the mass-in-box interface of each component, or null if not offered
    vector<char> _dustv;                // true for a dust component, which uses mass rather than number
    vector<int> _noIndices;             // empty list of indices returned when there is no medium system

    // data members describing the current node
    Box _box;
    int _level{0};

    // cached information for the current node; the flags are cleared by reset()
    bool _hasPositions{false};
    vector<Position> _positionv;                                  // the sample positions
    vector<char> _hasSamplesv;                                    // indexed on h
    vector<vector<double>> _samplesv;                             // indexed on h
    vector<char> _hasAmountv;                                     // indexed on h
    vector<double> _amountv;                                      // indexed on h
    std::array<bool, 3> _hasTypeSamplesv{{false, false, false}};  // indexed on material type
    std::array<vector<double>, 3> _typeSamplesv;                  // indexed on material type
};

//////////////////////////////////////////////////////////////////////

#endif
