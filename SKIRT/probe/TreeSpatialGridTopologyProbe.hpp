/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TREESPATIALGRIDTOPOLOGYPROBE_HPP
#define TREESPATIALGRIDTOPOLOGYPROBE_HPP

#include "SpecialtyProbe.hpp"

////////////////////////////////////////////////////////////////////

/** TreeSpatialGridTopologyProbe outputs a text data file representing the topology (i.e. the
    subdivision structure) of the tree grid used to discretize the spatial domain in the
    simulation. The resulting data file can be loaded by a TopologyTreePolicy instance to recreate
    a grid in another simulation that is identical to the one in the current simulation, or to
    refine that grid further with additional subdivision policies. This can be useful in situations
    where multiple simulations are being performed on input models with the same (or a very
    similar) spatial distribution of the transfer medium. Constructing a tree grid based on the
    medium density distribution takes time (because of the large number of density samples
    required), and the tree is likely to differ slightly between various runs (because the density
    is sampled at random positions, and the random number sequence varies because of
    parallelization). Loading the tree from the topology data file is much faster and guarantees
    that all simulations use identical tree grids. Note that the topology data stored in the file
    is scale-free; the simulation loading the topology must define the extent of the spatial
    domain.

    The output file is called <tt>prefix_probe_treetop.dat</tt>. It is a column text file with a
    single dimensionless column (see the ColumnOutFile class). After a brief descriptive header,
    which includes a column information line, the file contains lines with just a single integer
    number. The first line specifies the number of children for each nonleaf node (2 for a binary
    tree, 8 for an octtree, or 0 if the root node is not subdivided). The second line contains 1 if
    the root node is subdivided, or 0 if not. The following lines similarly contain 1 or 0
    indicating subdivision for any children of the preceding node, recursively, in a depth-first
    traversal of the tree. The children of a node are listed in the order defined by the
    BinTreeSpatialGrid and OctTreeSpatialGrid classes. Storing the number of children in the same
    column as the subdivision flags keeps the format of the files written by older versions, which
    lacked the column information line.

    If the output path of the simulation includes an HDF5 file, the topology is written as a text
    column file bundle in that file instead, with a single dataset holding the same values.

    By default, the probe is performed after setup, recording the topology of the grid constructed
    by the tree policies. If the grid is refined dynamically during the iterations of the
    simulation (see the DynamicRefinementOptions class), the \em probeAfter property can be set to
    record the topology of the refined grid at the end of the simulation run, or after each primary
    or secondary emission iteration. In the latter case, the name of the output file includes the
    iteration index, as in <tt>prefix_probe001_treetop.dat</tt>. A subsequent simulation can load
    the topology of the refined grid with the TopologyTreePolicy, and then starts from the refined
    grid. The reloaded grid has the same cells as the refined grid, but not in the same order,
    because the TopologyTreePolicy numbers the cells level by level, as for any tree, whereas
    dynamic refinement keeps the existing cell indices and appends the new cells. */
class TreeSpatialGridTopologyProbe : public SpecialtyProbe
{
    /** The enumeration type indicating when probing occurs. */
    ENUM_DEF(ProbeAfter, Setup, Run, Primary, Secondary)
        ENUM_VAL(ProbeAfter, Setup, "after setup")
        ENUM_VAL(ProbeAfter, Run, "after the complete simulation run")
        ENUM_VAL(ProbeAfter, Primary, "after each iteration over primary emission")
        ENUM_VAL(ProbeAfter, Secondary, "after each iteration over secondary emission")
    ENUM_END()

    ITEM_CONCRETE(TreeSpatialGridTopologyProbe, SpecialtyProbe,
                  "specialty: data file representing the topology of the tree spatial grid")
        ATTRIBUTE_TYPE_DISPLAYED_IF(TreeSpatialGridTopologyProbe, "Level3&TreeSpatialGrid")

        PROPERTY_ENUM(probeAfter, ProbeAfter, "perform the probe after")
        ATTRIBUTE_DEFAULT_VALUE(probeAfter, "Setup")
        ATTRIBUTE_DISPLAYED_IF(probeAfter, "IteratePrimary|IterateSecondary")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns an enumeration indicating when probing for this probe should be
        performed corresponding to the configured value of the \em probeAfter property. */
    When when() const override;

    //======================== Other Functions =======================

protected:
    /** This function performs probing. */
    void probe() override;
};

////////////////////////////////////////////////////////////////////

#endif
