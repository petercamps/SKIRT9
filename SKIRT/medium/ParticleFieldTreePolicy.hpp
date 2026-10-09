/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PARTICLEFIELDTREEPOLICY_HPP
#define PARTICLEFIELDTREEPOLICY_HPP

#include "SmoothingKernel.hpp"
#include "TreePolicy.hpp"
class ParticleSnapshot;

//////////////////////////////////////////////////////////////////////

/** ParticleFieldTreePolicy is a tree subdivision policy based on a scalar field defined by a set
    of smoothed particles imported from a text column file. A node is subdivided if it contains
    more than a given fraction of the total "mass" of the field. The field can reflect a property
    of the model other than the media actually used in the simulation, or the particles can be
    positioned "by hand" at strategic positions to obtain a higher resolution in the corresponding
    regions.

    The input file contains a particle on each line, with columns for the \f$x\f$, \f$y\f$, and
    \f$z\f$ coordinates of the particle's position, its smoothing length \f$h\f$, and its mass
    \f$M\f$. The default units are pc for the coordinates and the smoothing length, and Msun for
    the mass; this can be overridden by providing column information in the header of the file,
    as described for the TextInFile class. Because the policy uses only the fraction of the total
    mass contained in a node, the masses can be given in arbitrary units (as long as a column
    header, if present, specifies a mass unit). The smoothing kernel used to spread the mass of a
    particle over its smoothing sphere is configurable; the default is a cubic spline kernel.

    The mass contained in a node is calculated by integrating the smoothing kernel of each
    overlapping particle over the volume of the node, if the smoothing kernel offers this
    capability. Otherwise, it is estimated from the density of the field sampled at the random
    positions shared by the policies evaluating the node (see the TreeNodeEvaluation class). */
class ParticleFieldTreePolicy : public TreePolicy
{
    ITEM_CONCRETE(ParticleFieldTreePolicy, TreePolicy,
                  "a tree subdivision policy using a field of smoothed particles loaded from file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ParticleFieldTreePolicy, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the smoothed particles")

        PROPERTY_DOUBLE(maxFraction, "the maximum fraction of the field's mass contained in each cell")
        ATTRIBUTE_MIN_VALUE(maxFraction, "]0")
        ATTRIBUTE_MAX_VALUE(maxFraction, "1e-2]")
        ATTRIBUTE_DEFAULT_VALUE(maxFraction, "1e-6")

        PROPERTY_ITEM(smoothingKernel, SmoothingKernel, "the kernel for interpolating the smoothed particles")
        ATTRIBUTE_DEFAULT_VALUE(smoothingKernel, "CubicSplineSmoothingKernel")
        ATTRIBUTE_DISPLAYED_IF(smoothingKernel, "Level2")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor deletes the particle snapshot. */
    ~ParticleFieldTreePolicy();

protected:
    /** This function imports the smoothed particles. It is performed after the smoothing kernel
        has been set up. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns true if the mass of the field contained in the specified node
        exceeds the configured fraction of the total mass, and false otherwise. */
    bool needsSubdivide(TreeNodeEvaluation& node) const override;

    //======================== Data Members ========================

private:
    // data members initialized by setupSelfAfter()
    ParticleSnapshot* _snapshot{nullptr};  // the imported particles
    bool _hasMassInBox{false};             // true if the mass in a node can be calculated exactly
    double _maxMass{0.};                   // the maximum mass of the field in a node
};

//////////////////////////////////////////////////////////////////////

#endif
