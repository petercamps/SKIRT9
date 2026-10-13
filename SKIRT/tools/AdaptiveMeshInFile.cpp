/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AdaptiveMeshInFile.hpp"
#include "FatalError.hpp"
#include "H5BundleR.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

AdaptiveMeshInFile::AdaptiveMeshInFile(const SimulationItem* item, const string& filename, const string& description)
    : ColumnInFile(item, filename, description)
{}

////////////////////////////////////////////////////////////////////

bool AdaptiveMeshInFile::readNonLeaf(int& nx, int& ny, int& nz)
{
    open();
    if (TextInFile* text = textFile()) return text->readNonLeaf(nx, ny, nz);

    // the position of the next node follows from the number of nonleaf nodes and leaf nodes read so far
    size_t node = _nextNonLeaf + numRowsRead();
    if (node >= _isLeaf.size() || _isLeaf[node]) return false;

    nx = _topology[3 * _nextNonLeaf];
    ny = _topology[3 * _nextNonLeaf + 1];
    nz = _topology[3 * _nextNonLeaf + 2];
    _nextNonLeaf++;
    return true;
}

////////////////////////////////////////////////////////////////////

void AdaptiveMeshInFile::openedBundle(const H5BundleR& bundle)
{
    // read the flags that distinguish leaf nodes from nonleaf nodes
    H5DatasetR isLeaf = bundle.getDataset("is_leaf");
    if (isLeaf.shape().size() != 1)
        throw FATALERROR("Dataset 'is_leaf' in bundle " + bundlePath() + " is not one-dimensional");
    _isLeaf = isLeaf.readBool();
    size_t numLeaves = std::count_if(_isLeaf.begin(), _isLeaf.end(), [](char flag) { return flag != 0; });
    size_t numNonLeaves = _isLeaf.size() - numLeaves;

    // read the subdivision counts for the nonleaf nodes
    H5DatasetR topology = bundle.getDataset("topology");
    vector<size_t> shape = topology.shape();
    if (shape.size() != 2 || shape[1] != 3)
        throw FATALERROR("Dataset 'topology' in bundle " + bundlePath() + " does not have three columns");
    if (shape[0] != numNonLeaves)
        throw FATALERROR("Dataset 'topology' in bundle " + bundlePath() + " has " + std::to_string(shape[0])
                         + " rows rather than one for each of the " + std::to_string(numNonLeaves) + " nonleaf nodes");
    _topology = topology.readInt();

    // verify that the nodes form a single complete tree: each node fills one open position, and a nonleaf node
    // opens a position for each of its children; the tree is complete when no open positions remain
    std::int64_t numOpen = 1;
    size_t nonLeaf = 0;
    for (char flag : _isLeaf)
    {
        if (numOpen == 0)
            throw FATALERROR("Adaptive mesh in bundle " + bundlePath()
                             + " has superfluous nodes after all nodes of the tree were specified");
        numOpen--;
        if (!flag)
        {
            std::int64_t nx = _topology[3 * nonLeaf];
            std::int64_t ny = _topology[3 * nonLeaf + 1];
            std::int64_t nz = _topology[3 * nonLeaf + 2];
            if (nx < 1 || ny < 1 || nz < 1)
                throw FATALERROR("Dataset 'topology' in bundle " + bundlePath()
                                 + " has a subdivision count that is not positive");
            numOpen += nx * ny * nz;
            nonLeaf++;
        }
    }
    if (numOpen != 0)
        throw FATALERROR("Adaptive mesh in bundle " + bundlePath()
                         + " ends before all nodes of the tree were specified");

    // each property of the leaf nodes has a row for each leaf node
    setNumRows(numLeaves);
}

////////////////////////////////////////////////////////////////////
