/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BINTREESPATIALGRID_HPP
#define BINTREESPATIALGRID_HPP

#include "BoxCellDensityMixIn.hpp"
#include "BoxSpatialGrid.hpp"

//////////////////////////////////////////////////////////////////////

/** BinTreeSpatialGrid is a concrete subclass of the BoxSpatialGrid class that represents a
    three-dimensional spatial grid with cuboidal cells organized in a binary tree, loaded from a
    topology data file. The tree's root node encloses the complete spatial domain. Each node is
    either a leaf node (an actual spatial cell), or is divided into two equal halves along a plane
    perpendicular to one of the coordinate axes, alternating between the x, y and z axes when
    descending the tree. The depth of the tree can vary from place to place.

    The topology data file has the same format as for the FileTreeSpatialGrid class. In the most
    common use case, it has been created in a previous simulation through the
    TreeSpatialGridTopologyProbe. After a brief descriptive header, the file contains lines with
    just a single integer number. The first line specifies the number of children for each nonleaf
    node (2 for a binary tree, or 0 if the root node is not subdivided; an octtree with 8 children
    is not supported by this class). The second line contains 1 if the root node is subdivided, or
    0 if not. The following lines similarly contain 1 or 0 indicating subdivision for any children
    of the preceding node, recursively, in a depth-first traversal of the tree. The topology data
    is scale-free; the extent of the spatial domain must be configured in the simulation loading
    the topology.

    Given the same topology file and spatial domain, this class constructs exactly the same tree
    grid as FileTreeSpatialGrid, with the same cell indices and identical cell boundaries, so that
    the two can be used interchangeably. The difference is purely internal: this class is
    optimized for speed and memory use rather than for generality. It does not derive from
    TreeSpatialGrid and its TreeNode hierarchy, but stores the tree in a single contiguous array of
    small nodes that refer to each other by index rather than by pointer, and that are handled
    without any virtual function calls. Each node holds its own extent, the index of its children,
    the cell index (for a leaf node), and the index of its neighbor across each of its six walls.
    The neighbor across a wall is the node at the same level in the tree that touches the wall if
    there is one; otherwise it is the leaf node (at a lower level) that covers the other side of
    the wall. Neighbor links are established for all nodes in a single pass through the tree after
    it has been loaded, rather than through incremental bookkeeping while nodes are subdivided.

    The path segment generator uses these links to move from a leaf node to the next one after
    crossing a wall. If the neighbor turns out to have children of its own (because the tree is
    more refined on the other side of the wall), the generator descends into the neighbor until it
    reaches the leaf containing the new position. Compared to searching a list of neighbors, this
    requires just a single memory access in the most common case, and it never depends on the
    number of cells adjacent to a wall. A top-down search starting at the root node is used to
    determine the initial cell and as a fall-back in rare cases where numerical inaccuracies would
    otherwise result in an inconsistent state, just like the implementation for other tree grids.

    The nodes of a large tree are usually not in the processor cache, so that waiting for the next
    node to arrive from memory dominates the cost of each step. To hide this latency, the
    generator asks the processor to prefetch the (up to) three nodes that the path can move to next
    while it is still working on the current node. This relies on a compiler-specific builtin
    function that is used only for compilers offering it (GCC, Clang and the Intel compilers); for
    other compilers prefetching is skipped.

    For paths that do not run exactly along the boundary between two cells, this class visits the
    same cells as FileTreeSpatialGrid, and the lengths of the path segments agree to within
    floating point rounding errors (this class multiplies with the reciprocal of the direction,
    which is calculated just once for each path, rather than dividing in each step). A path that
    does run exactly along a cell boundary (for example, a path parallel to a coordinate axis
    through a position on a splitting plane) is ambiguous: the cell before and the cell after the
    boundary are equally valid choices. In such cases the two classes may select different,
    equally valid cells. This class consistently selects the cell on the upper side of a splitting
    plane, as it does when locating the cell that contains a given position.

    Because it is not a TreeSpatialGrid, this class cannot be used with the
    TreeSpatialGridTopologyProbe. */
class BinTreeSpatialGrid : public BoxSpatialGrid, public BoxCellDensityMixIn
{
    ITEM_CONCRETE(BinTreeSpatialGrid, BoxSpatialGrid, "a binary tree spatial grid loaded from a topology data file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(BinTreeSpatialGrid, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the tree topology data")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function reads the topology of the tree from the configured data file, constructs the
        internal representation of the tree (including the links between neighboring nodes), and
        logs some details on the number of cells in the tree. */
    void setupSelfAfter() override;

    //======================== Other Functions =======================

public:
    /** This function returns the number of cells in the grid. */
    int numCells() const override;

    /** This function returns the box defining the cell with index \f$m\f$, as required by the
        BoxCellDensityMixIn class. */
    Box cellBox(int m) const override;

    /** This function returns the volume of the cell with index \f$m\f$, calculated as \f$V =
        \Delta x\, \Delta y\, \Delta z\f$. */
    double volume(int m) const override;

    /** This function returns the actual diagonal of the cell with index \f$m\f$, calculated as
        \f$d = \sqrt{ (\Delta x)^2 + (\Delta y)^2 + (\Delta z)^2 }\f$. */
    double diagonal(int m) const override;

    /** This function returns the index of the cell that contains the position \f${\bf{r}}\f$. The
        search algorithm starts at the root node and selects the child node that contains the
        position. This procedure is repeated until the node is childless, i.e. until it is a leaf
        node that corresponds to an actual spatial cell. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$, calculated as
        \f[ \begin{split} x &= x_{\text{min}} + \frac12\, \Delta x \\ y &= y_{\text{min}} +
        \frac12\, \Delta y \\ z &= z_{\text{min}} + \frac12\, \Delta z \end{split} \f] */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$, calculated as
        \f[ \begin{split} x &= x_{\text{min}} + {\cal{X}}_1\, \Delta x \\ y &= y_{\text{min}} +
        {\cal{X}}_2\, \Delta y \\ z &= z_{\text{min}} + {\cal{X}}_3\, \Delta z \end{split} \f] with
        \f${\cal{X}}_1\f$, \f${\cal{X}}_2\f$ and \f${\cal{X}}_3\f$ three uniform deviates. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for this grid, implemented as a private
        PathSegmentGenerator subclass. The algorithm is the same as the one used for other tree
        grids. It determines the cell that contains the starting position, and calculates the first
        wall of the cell that will be crossed. The pathlength \f$\Delta s\f$ is determined and the
        current position is moved to a new position along this path, a tiny fraction further than
        \f$\Delta s\f$, \f[ \begin{split} x_{\text{new}} &= x_{\text{current}} + (\Delta s +
        \epsilon)\,k_x \\ y_{\text{new}} &= y_{\text{current}} + (\Delta s + \epsilon)\,k_y \\
        z_{\text{new}} &= z_{\text{current}} + (\Delta s + \epsilon)\,k_z \end{split} \f] where \f[
        \epsilon = 10^{-12} \sqrt{x_{\text{max}}^2 + y_{\text{max}}^2 + z_{\text{max}}^2} \f] By
        adding this small extra bit, we ensure that the new position is now within the next cell,
        and we can repeat this exercise. This loop is terminated when the next position is outside
        the grid.

        The difference with the other tree grids is in the way the next cell is determined: this
        implementation follows the link to the neighbor across the crossed wall, and if needed
        descends from there to the leaf node that contains the new position, as described in the
        class header. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

protected:
    /** This function writes the intersection of the grid with the xy plane to the specified
        SpatialGridPlotFile object. */
    void write_xy(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid with the xz plane to the specified
        SpatialGridPlotFile object. */
    void write_xz(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid with the yz plane to the specified
        SpatialGridPlotFile object. */
    void write_yz(SpatialGridPlotFile* outfile) const override;

    /** This function writes 3D information for the cells up to a certain level in the grid
        structure to the specified SpatialGridPlotFile object. The output is limited to some
        predefined number of cells to keep the line density in the output plot within reason. */
    void write_xyz(SpatialGridPlotFile* outfile) const override;

    /** This function is used by the interface() function to ensure that the receiving item can
        actually offer the specified interface. If the requested interface is the
        DensityInCellInterface, the implementation in this class returns the value returned by the
        BoxCellDensityMixIn::offersInterface() function. For other requested interfaces, the
        function invokes its counterpart in the base class. */
    bool offersInterface(const std::type_info& interfaceTypeInfo) const override;

    //======================== Data Members ========================

private:
    // the implementation of the tree data structure is hidden in the source file; it is held by a shared
    // pointer (rather than a unique pointer) so that this class can be constructed without a complete definition
    class Tree;

    // data members initialized during setup
    double _eps{0.};                    // a small fraction relative to the spatial extent of the grid
    std::shared_ptr<const Tree> _tree;  // the tree data structure

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

//////////////////////////////////////////////////////////////////////

#endif
