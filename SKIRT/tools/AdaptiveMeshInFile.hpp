/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ADAPTIVEMESHINFILE_HPP
#define ADAPTIVEMESHINFILE_HPP

#include "ColumnInFile.hpp"

////////////////////////////////////////////////////////////////////

/** This class allows reading the description of an adaptive mesh from a user input file. An
    adaptive mesh is a hierarchical tree of nodes, in which each nonleaf node is subdivided into a
    regular grid of child nodes, and each leaf node (a cell) carries a number of properties. The
    mesh can be provided as a regular column text file, or as an AMR text file bundle in an HDF5
    file. The class inherits the ColumnInFile class, which determines where the input is read from
    and which provides the functions to declare and read the properties of the leaf nodes, and it
    adds the readNonLeaf() function to read the structure of the tree.

    In both formats, the nodes of the tree are listed in the same order: the root node first, and
    the children of a nonleaf node immediately after that node, recursively, with the first index
    of the grid of children varying most rapidly. A client reads the nodes in that order: it calls
    the readNonLeaf() function, and if that function returns false because the next node is a
    leaf node, it calls one of the readRow() functions to obtain the properties of the leaf node.

    Column text file
    ----------------

    In a column text file, a nonleaf node is represented by a line that starts with an exclamation
    mark, followed by the number of subdivisions in each of the three spatial directions. A leaf
    node is represented by a regular line with a value for each property. See the TextInFile class
    for more information.

    AMR text file bundle
    --------------------

    An AMR text file bundle holds the following datasets:

    - \em is_leaf: a one-dimensional dataset of Booleans with a value for each node of the tree, in
      the order described above: true for a leaf node and false for a nonleaf node.

    - \em topology: a two-dimensional dataset of integers with a row for each nonleaf node, in the
      order in which these nodes occur, and three columns that specify the number of subdivisions
      of the node in the x, y, and z direction.

    - a one-dimensional dataset for each property of the leaf nodes, with a value for each leaf
      node in the order in which these nodes occur, and with the attributes described for a text
      column file bundle in the ColumnInFile class.

    When the bundle is opened, this class verifies that the nodes form a single, complete tree. */
class AdaptiveMeshInFile final : public ColumnInFile
{
    //=============== Construction - Destruction  ==================

public:
    /** The constructor records its arguments without opening anything, as described for the
        ColumnInFile constructor. */
    AdaptiveMeshInFile(const SimulationItem* item, const string& filename, const string& description);

    //========================== Reading ===========================

public:
    /** This function attempts to read a nonleaf node. If the next node in the file is a nonleaf
        node, the function stores the number of subdivisions of that node in each spatial
        direction in its arguments, advances to the next node, and returns true. Otherwise, i.e. if
        the next node is a leaf node or if there are no more nodes, the function returns false
        without advancing, and the contents of its arguments is undefined. In the first case, the
        client can then call one of the readRow() functions to read the properties of the leaf
        node. This function has the same signature and behavior as the TextInFile::readNonLeaf()
        function. */
    bool readNonLeaf(int& nx, int& ny, int& nz);

    //======================== Private helpers ========================

protected:
    /** This function reads the datasets that describe the structure of the tree from the specified
        bundle, verifies them, and informs the base class of the number of leaf nodes. */
    void openedBundle(const H5BundleR& bundle) override;

    //======================== Data Members ========================

private:
    // the structure of the tree, if the input is being read from a bundle
    vector<char> _isLeaf;    // for each node, nonzero for a leaf node and zero for a nonleaf node
    vector<int> _topology;   // for each nonleaf node, the number of subdivisions in the x, y, and z direction
    size_t _nextNonLeaf{0};  // the zero-based index of the next nonleaf node to be read
};

////////////////////////////////////////////////////////////////////

#endif
