/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TREESPATIALGRID_HPP
#define TREESPATIALGRID_HPP

#include "BoxCellDensityMixIn.hpp"
#include "BoxSpatialGrid.hpp"
#include "TreePolicy.hpp"
#include <array>
#include <deque>
class TextOutFile;

//////////////////////////////////////////////////////////////////////

/** TreeSpatialGrid is an abstract subclass of the BoxSpatialGrid class, and represents
    three-dimensional spatial grids with cuboidal cells organized in a hierarchical tree. The
    tree's root node encloses the complete spatial domain, and nodes on subsequent levels
    recursively divide space into ever finer nodes. The depth of the tree can vary from place to
    place. The leaf nodes (those that are not further subdivided) are the actual spatial cells.

    <b>Construction</b>

    The tree is constructed according to a list of user-configured subdivision policies (see the
    TreePolicy class) within the configured range of subdivision levels: nodes are always
    subdivided up to the minimum level, and never beyond the maximum level. In between, a node is
    subdivided as soon as one of the policies asks for it, so that the policies can be combined
    freely. The policies are evaluated in the order in which they are listed, and evaluation stops
    as soon as a policy asks for subdivision. If the list is empty, the tree is subdivided to the
    minimum level everywhere.

    The tree is constructed level by level. Starting from a list holding just the root node, each
    pass evaluates every node at the current level, and subdivides the nodes that need it,
    appending their children to the end of the same list. The newly appended range of nodes thus
    forms the next level. Because subdivision only ever appends, the position of a node in the list
    serves as a stable, level-ordered identifier. Evaluating whether a node needs subdivision is
    read-only and can be expensive (for example, it may require sampling the density of the media),
    so it is performed in parallel for all nodes at a level. The policies evaluating a node share
    the properties of the media in the node, which are calculated only once (see the
    TreeNodeEvaluation class). Subdividing the flagged nodes is performed sequentially.

    The way in which a node is subdivided depends on the type of tree, and is implemented by a
    subclass: an octtree (8 children per node) or a binary tree (2 children per node). In both
    cases, the children of a node are consecutive in the list.

    <b>Representation</b>

    Once construction is complete, the tree is stored as a flat array of small, fixed-size nodes
    that refer to each other by index into the array rather than by pointer. Each node holds its own
    extent, the index of its first child (the other children follow it directly), its cell index
    (for a leaf node), and the index of its neighbor across each of its six walls. The neighbor
    across a wall is the node at the same level in the tree that touches the wall, if there is one;
    otherwise it is the leaf node (at a lower level) that covers the other side of the wall. The
    neighbor links are established by the subclass in a single top-down pass over the complete
    tree. The cells are the leaf nodes, numbered in the order in which they occur in the array.

    This base class implements the functions that depend only on the extent of the cells, and the
    output of the grid structure and topology. The subclass implements the functions that depend on
    the tree type, including the location of the cell containing a given position and the
    generation of path segments, which follows the neighbor links. */
class TreeSpatialGrid : public BoxSpatialGrid, public BoxCellDensityMixIn
{
    ITEM_ABSTRACT(TreeSpatialGrid, BoxSpatialGrid, "a hierarchical tree spatial grid")

        PROPERTY_ITEM_LIST(policies, TreePolicy, "the tree subdivision policies")
        ATTRIBUTE_DEFAULT_VALUE(policies, "DensityTreePolicy")
        ATTRIBUTE_REQUIRED_IF(policies, "false")

        PROPERTY_INT(minLevel, "the minimum level of grid refinement")
        ATTRIBUTE_MIN_VALUE(minLevel, "0")
        ATTRIBUTE_MAX_VALUE(minLevel, "99")
        ATTRIBUTE_DEFAULT_VALUE(minLevel, "OctTreeSpatialGrid:3;BinTreeSpatialGrid:9")

        PROPERTY_INT(maxLevel, "the maximum level of grid refinement")
        ATTRIBUTE_MIN_VALUE(maxLevel, "0")
        ATTRIBUTE_MAX_VALUE(maxLevel, "99")
        ATTRIBUTE_DEFAULT_VALUE(maxLevel, "OctTreeSpatialGrid:7;BinTreeSpatialGrid:21")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function verifies that the maximum level is not below the minimum level. */
    void setupSelfBefore() override;

    /** This function constructs the tree as described in the class header, converts it to the
        flat array representation, asks the subclass to establish the neighbor links, determines
        the cell indices, and logs some details on the number of cells in the tree. */
    void setupSelfAfter() override;

    //======================== Nodes =======================

public:
    /** The Node class represents a node in the tree, which is either a leaf node (corresponding to
        a spatial cell) or a nonleaf node with children. The nodes are stored by value in a single
        contiguous array, and they refer to each other by index into this array. The children of a
        nonleaf node are consecutive in the array, so a node stores just the index of its first
        child. The class has no virtual functions.

        The class is public only so that the implementation files of the subclasses can use it in
        local functions. Other classes cannot obtain the nodes of a tree, because the nodes()
        function is protected.

        The six walls of a node are numbered such that wall \f$2a\f$ is the lower and wall
        \f$2a+1\f$ the upper wall perpendicular to axis \f$a\f$, with \f$a=0,1,2\f$ for x, y, and
        z. A node stores the coordinate of each wall along the corresponding axis, and the index of
        its neighbor across each wall (or -1 for a wall on the boundary of the domain). It also
        stores the axis perpendicular to the splitting plane, which is used by binary trees only.
        Indices are 32-bit integers to keep the node small: its size is 88 bytes. */
    class Node
    {
    public:
        /** This constructor creates a node with the specified extent and level, without any
            children or neighbors. The splitting axis is set to the level modulo three. */
        Node(const Box& extent, int level)
            : _wallv{{extent.xmin(), extent.xmax(), extent.ymin(), extent.ymax(), extent.zmin(), extent.zmax()}},
              _level(level), _axis(level % 3)
        {}

        /** This function returns the extent of the node as a box. */
        Box extent() const { return Box(_wallv[0], _wallv[2], _wallv[4], _wallv[1], _wallv[3], _wallv[5]); }

        /** This function returns the level of the node in the tree, with level zero for the root
            node. */
        int level() const { return _level; }

        /** This function returns the axis (0=x, 1=y, 2=z) perpendicular to the plane along which a
            nonleaf node in a binary tree is split. */
        int axis() const { return _axis; }

        /** This function returns the coordinate of the specified wall (0-5) along the axis
            perpendicular to that wall. */
        double wall(int w) const { return _wallv[w]; }

        /** This function returns the coordinate of the center of the node along the specified axis
            (0-2), i.e. the coordinate of a splitting plane perpendicular to that axis. */
        double center(int a) const { return 0.5 * (_wallv[2 * a] + _wallv[2 * a + 1]); }

        /** This function returns true if the specified position is inside the node, borders
            included. */
        bool contains(double x, double y, double z) const
        {
            return x >= _wallv[0] && x <= _wallv[1] && y >= _wallv[2] && y <= _wallv[3] && z >= _wallv[4]
                   && z <= _wallv[5];
        }

        /** This function returns true if the top-down search for the leaf node containing the
            specified position, which assigns a position on a splitting plane to the child on the
            upper side, would end in this leaf node. In other words, it returns true if the
            position is inside the node, including its lower walls and excluding its upper walls,
            except for upper walls on the boundary of the domain (i.e. without a neighbor), which
            are included. */
        bool owns(double x, double y, double z) const
        {
            return x >= _wallv[0] && (x < _wallv[1] || (x == _wallv[1] && _neighborv[1] < 0)) && y >= _wallv[2]
                   && (y < _wallv[3] || (y == _wallv[3] && _neighborv[3] < 0)) && z >= _wallv[4]
                   && (z < _wallv[5] || (z == _wallv[5] && _neighborv[5] < 0));
        }

        /** This function returns true if the node is a leaf node, i.e. a spatial cell without
            children. */
        bool isLeaf() const { return _child < 0; }

        /** This function returns the index of the first child of a nonleaf node; the other
            children follow it directly. */
        int child() const { return _child; }

        /** This function returns the cell index for a leaf node, or -1 for a nonleaf node. */
        int cell() const { return _cell; }

        /** This function returns the index of the neighbor across the specified wall (0-5), or -1
            if there is none. */
        int neighbor(int w) const { return _neighborv[w]; }

        /** This function sets the index of the first child. */
        void setChild(int child) { _child = child; }

        /** This function sets the cell index. */
        void setCell(int cell) { _cell = cell; }

        /** This function sets the index of the neighbor across the specified wall (0-5). */
        void setNeighbor(int w, int neighbor) { _neighborv[w] = neighbor; }

        /** This function sets the coordinate of the specified wall (0-5). */
        void setWall(int w, double coordinate) { _wallv[w] = coordinate; }

    private:
        std::array<double, 6> _wallv;                             // xmin, xmax, ymin, ymax, zmin, zmax
        std::array<int, 6> _neighborv{{-1, -1, -1, -1, -1, -1}};  // neighbor across each wall
        int _child{-1};                                           // index of the first child; -1 for a leaf
        int _cell{-1};                                            // cell index; -1 for a nonleaf node
        int _level{0};                                            // level of the node in the tree
        int _axis{0};                                             // splitting axis for a binary tree
    };

    //=========== Functions implemented by subclasses ===========

public:
    /** This function must be implemented in a subclass to return the number of children of a
        nonleaf node: 2 for a binary tree and 8 for an octtree. */
    virtual int numChildren() const = 0;

    /** This function must be implemented in a subclass to return the extent of the child with
        index \f$c\f$ (from zero to numChildren()-1) of a node with the specified extent and level
        in the tree. The tree uses this function to subdivide its nodes. Together with the
        childIndex() function, it also allows a subdivision policy to mirror the structure of the
        tree, for example to descend a recorded tree topology. */
    virtual Box childExtent(const Box& extent, int level, int c) const = 0;

    /** This function must be implemented in a subclass to return the index (from zero to
        numChildren()-1) of the child that contains the specified position, for a node with the
        specified extent and level in the tree. A position on a splitting plane is assigned to the
        child on the upper side of the plane. The position is assumed to be inside the node. */
    virtual int childIndex(const Box& extent, int level, Vec position) const = 0;

protected:
    /** This function must be implemented in a subclass to establish the neighbor links for all
        nodes in the specified array, as described in the class header. On entry, the nodes have
        their extent, level, and children, but no neighbors. The nodes are ordered such that each
        node comes after its parent, and the root node has index zero. */
    virtual void linkNeighbors(vector<Node>& nodes) const = 0;

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

    /** This function returns the central location of the cell with index \f$m\f$, calculated as
        \f[ \begin{split} x &= x_{\text{min}} + \frac12\, \Delta x \\ y &= y_{\text{min}} +
        \frac12\, \Delta y \\ z &= z_{\text{min}} + \frac12\, \Delta z \end{split} \f] */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$, calculated as
        \f[ \begin{split} x &= x_{\text{min}} + {\cal{X}}_1\, \Delta x \\ y &= y_{\text{min}} +
        {\cal{X}}_2\, \Delta y \\ z &= z_{\text{min}} + {\cal{X}}_3\, \Delta z \end{split} \f] with
        \f${\cal{X}}_1\f$, \f${\cal{X}}_2\f$ and \f${\cal{X}}_3\f$ three uniform deviates. */
    Position randomPositionInCell(int m) const override;

    /** This function writes the topology of the tree to the specified text file in a simple,
        proprietary format. After a brief descriptive header, it writes lines that each contain
        just a single integer number. The first line specifies the number of children for each
        nonleaf node (2 for a binary tree, 8 for an octtree, or 0 if the root node has not been
        subdivided). The second line contains 1 if the root node is subdivided, or 0 if not. The
        following lines similarly contain 1 or 0 indicating subdivision for any children of the
        preceding node, recursively, in a depth-first traversal of the tree. */
    void writeTopology(TextOutFile* outfile) const;

protected:
    /** This function returns the array of nodes; the first node is the root node. It is intended
        for use by subclasses, after setup has been completed. */
    const vector<Node>& nodes() const { return _nodev; }

    /** This function returns, for each cell index, the index of the corresponding leaf node in the
        array returned by the nodes() function. It is intended for use by subclasses, after setup
        has been completed. */
    const vector<int>& cellNodeIndices() const { return _idv; }

    /** This function returns a small distance relative to the spatial extent of the grid, used by
        the path segment generator to move a position just across a cell wall. */
    double eps() const { return _eps; }

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

private:
    /** This function constructs the tree level by level according to the configured policies and
        subdivision levels, and returns the resulting list of nodes. */
    std::deque<Node> constructTree() const;

    /** This function returns a reference to the leaf node corresponding to the cell with index
        \f$m\f$. */
    const Node& cellNode(int m) const { return _nodev[_idv[m]]; }

    //======================== Data Members ========================

private:
    // data members initialized during setup
    double _eps{0.};      // a small fraction relative to the spatial extent of the grid
    vector<Node> _nodev;  // all nodes in the tree; the first one is the root node
    vector<int> _idv;     // index in _nodev of the leaf node for each cell (i.e. for each cell index)
};

//////////////////////////////////////////////////////////////////////

#endif
