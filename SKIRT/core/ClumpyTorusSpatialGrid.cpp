/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ClumpyTorusSpatialGrid.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "PathSegmentGenerator.hpp"
#include "Random.hpp"
#include "SpatialGridPlotFile.hpp"
#include "TextInFile.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // This function returns the square of a value
    double square(double v)
    {
        return v * v;
    }

    // This function returns the cube of a value
    double cube(double v)
    {
        return v * v * v;
    }

    // This function returns the volume of a sphere with given radius
    double volumeSphere(double radius)
    {
        return (4. / 3.) * M_PI * cube(radius);
    }

    // This function returns the volume of a torus with given opening angle and radii
    double volumeTorus(double delta, double r1, double r2)
    {
        return (4. / 3.) * M_PI * std::sin(delta) * (cube(r2) - cube(r1));
    }

    // This function returns true if the sphere with given center and radius is fully inside the torus
    // with given opening angle and radii (no straggling nor touching allowed)
    bool isSphereInTorus(double x, double y, double z, double r, double delta, double r1, double r2)
    {
        // sphere must be strictly outside the inner radius and inside the outer radius
        const double rho = std::sqrt(square(x) + square(y) + square(z));
        if (rho - r <= r1 || rho + r >= r2) return false;

        // sphere must be strictly inside the conical wedge |theta| < delta
        const double s = std::sqrt(square(x) + square(y));
        const double theta = std::atan2(z, s);    // in [-pi/2, pi/2]
        const double alpha = std::asin(r / rho);  // half-angle subtended, in (0, pi/2)
        if (std::fabs(theta) + alpha >= delta) return false;

        return true;
    }

    // This function returns true if the given position is inside the torus
    // with given opening angle and radii (no touching allowed)
    bool isPositionInTorus(Position p, double delta, double r1, double r2)
    {
        // position must be strictly outside the inner radius and inside the outer radius
        const double rho = p.radius();
        if (rho <= r1 || rho >= r2) return false;

        // position must be strictly inside the conical wedge |theta| < delta
        const double s = p.cylRadius();
        const double theta = std::atan2(p.height(), s);  // in [-pi/2, pi/2]
        if (std::fabs(theta) >= delta) return false;

        return true;
    }

    // This function returns true if the specified spheres overlap or touch, and false otherwise
    double doSpheresOverlap(Position center1, double radius1, Position center2, double radius2)
    {
        return (center1 - center2).norm2() <= square(radius1 + radius2);
    }

    // This function returns a random position in a uniform sphere with given center and radius
    Position randomPositionInSphere(Random* random, Position center, double radius)
    {
        double r = std::cbrt(random->uniform());
        Direction k = random->direction();
        return Position(center + k * r * radius);
    }

    // This function returns a random position in a uniform torus with given opening angle and radii
    Position randomPositionInTorus(Random* random, double delta, double r1, double r2)
    {
        // radius
        double X = random->uniform();
        double r = std::cbrt((1. - X) * cube(r1) + X * cube(r2));

        // inclination
        X = random->uniform();
        double costheta = (1. - 2. * X) * std::sin(delta);
        double theta = std::acos(costheta);

        // azimuth
        X = random->uniform();
        double phi = 2. * M_PI * X;

        return Position(r, theta, phi, Position::CoordinateSystem::SPHERICAL);
    }
}

////////////////////////////////////////////////////////////////////

// data type defining a custom bounding volume hierarchy for clumps
class ClumpyTorusSpatialGrid::BVH
{
private:
    class Node
    {
    public:
        Box box;         // union of the bounding boxes of all entities in this node's subtree
        int left;        // index into _nodes of the left child, or -1 if this is a leaf
        int right;       // index into _nodes of the right child, or -1 if this is a leaf
        int firstIndex;  // for a leaf: start of this leaf's range in _index
        int numIndices;  // for a leaf: number of entities in this leaf; 0 for an interior node

        Node(const Box& box_, int left_, int right_, int firstIndex_, int numIndices_)
            : box(box_), left(left_), right(right_), firstIndex(firstIndex_), numIndices(numIndices_)
        {}

        bool isLeaf() const { return numIndices > 0; }
    };
    // flat array of BVH nodes; index 0 is the root, meaningful only when _numEntities > 0
    vector<Node> _nodes;

    // entity indices, reordered during the build so that the entities of any single leaf
    // occupy the contiguous range [firstIndex, firstIndex+numIndices) of this array
    vector<int> _index;

    // pointer to array with clumps
    const vector<Clump>* _clumps;

    // private recursive build function
    int buildRecursive(int begin, int end) { /*TO DO*/ }

public:
    BVH() {}
    void loadClumps(const vector<Clump>& clumps)  { /*TO DO*/ }
    vector<int> allDisjointClumps() const  { /*TO DO*/ }
    int anyClumpContaining(Vec bfr) const  { /*TO DO*/ }
};

////////////////////////////////////////////////////////////////////

void ClumpyTorusSpatialGrid::setupSelfBefore()
{
    SpatialGrid::setupSelfBefore();
    Log* log = find<Log>();

    // verify configuration values
    if (_maxRadius <= _minRadius) throw FATALERROR("Maximum radius must be larger than minimum radius");

    // read the file defining the spherical clumps, rejecting those outside the domain
    int numOutside = 0;
    {
        TextInFile infile(this, _filename, "clump centers and radii");
        infile.addColumn("position x", "length", "pc");
        infile.addColumn("position y", "length", "pc");
        infile.addColumn("position z", "length", "pc");
        infile.addColumn("radius r", "length", "pc");
        double x, y, z, r;
        while (infile.readRow(x, y, z, r))
        {
            if (isSphereInTorus(x, y, z, r, _openingAngle, _minRadius, _maxRadius))
            {
                _clumps.emplace_back(x, y, z, r);
            }
            else
                numOutside++;
        }
        infile.close();
    }

    // build a BVH over all clumps inside the torus, including mutually overlapping ones
    _bvh = new BVH;
    log->info("Constructing bounding volume hierarchy for " + std::to_string(_clumps.size()) + " clumps...");
    _bvh->loadClumps(_clumps);

    // check for overlapping clumps using the BVH
    log->info("Checking for overlapping clumps...");
    vector<int> disjointIndices = _bvh->allDisjointClumps();
    int numOverlapping = _clumps.size() - disjointIndices.size();

    // if there are any overlapping clumps, rebuild the list of clumps and rebuild the BVH
    if (numOverlapping > 0)
    {
        log->info("Rebuilding bounding volume hierarchy for non-overlapping clumps...");
        vector<Clump> originalClumps;
        originalClumps.swap(_clumps);
        for (int m : disjointIndices) _clumps.emplace_back(originalClumps[m]);
        _bvh->loadClumps(_clumps);
    }

    // remember the final number of clumps, which also serves as the cell index for the torus
    _numClumps = _clumps.size();

    // inform the user
    log->info("Summary:");
    log->info("  Clumps read from file:     " + std::to_string(_numClumps + numOutside + numOverlapping));
    log->info("  Rejected (outside domain): " + std::to_string(numOutside));
    log->info("  Rejected (overlapping):    " + std::to_string(numOverlapping));
    log->info("  Remaining in spatial grid: " + std::to_string(_numClumps));
}

//////////////////////////////////////////////////////////////////////

ClumpyTorusSpatialGrid::~ClumpyTorusSpatialGrid()
{
    delete _bvh;
}

//////////////////////////////////////////////////////////////////////

int ClumpyTorusSpatialGrid::dimension() const
{
    return 3;
}

//////////////////////////////////////////////////////////////////////

int ClumpyTorusSpatialGrid::numCells() const
{
    return _numClumps + 1;
}

//////////////////////////////////////////////////////////////////////

Box ClumpyTorusSpatialGrid::boundingBox() const
{
    const double Rmax = _maxRadius;
    const double zmax = _maxRadius * std::sin(_openingAngle);
    return Box(-Rmax, -Rmax, -zmax, Rmax, Rmax, zmax);
}

//////////////////////////////////////////////////////////////////////

double ClumpyTorusSpatialGrid::volume(int m) const
{
    if (m >= 0 && m < _numClumps) return volumeSphere(_clumps[m].radius());
    if (m == _numClumps)
    {
        double volumeClumps = 0.;
        for (const auto& clump : _clumps) volumeClumps += volumeSphere(clump.radius());
        return volumeTorus(_openingAngle, _minRadius, _maxRadius) - volumeClumps;
    }
    return 0.;
}

//////////////////////////////////////////////////////////////////////

double ClumpyTorusSpatialGrid::diagonal(int m) const
{
    if (m >= 0 && m < _numClumps) return 2. * _clumps[m].radius();
    if (m == _numClumps) return 2. * (_maxRadius - _minRadius);  // use a simplistic estimate
    return 0.;
}

//////////////////////////////////////////////////////////////////////

int ClumpyTorusSpatialGrid::cellIndex(Position bfr) const
{
    auto contains = [this](int m, Vec p) { return (p - _clumps[m].center()).norm2() <= square(_clumps[m].radius()); };
    int m = _bvh->anyClumpContaining(bfr);
    if (m >= 0) return m;

    // return the torus index if the position is in the torus but not in any of the clumps
    if (isPositionInTorus(bfr, _openingAngle, _minRadius, _maxRadius)) return _numClumps;
    return -1;
}

//////////////////////////////////////////////////////////////////////

Position ClumpyTorusSpatialGrid::centralPositionInCell(int m) const
{
    if (m >= 0 && m < _numClumps) return _clumps[m].center();
    if (m == _numClumps)
    {
        // try a position halfway the x-axis
        Position candidate(0.5 * (_minRadius + _maxRadius), 0., 0.);
        if (cellIndex(candidate) == _numClumps) return candidate;

        // if that position is inside a clump, return a random position instead
        return randomPositionInCell(_numClumps);
    }
    return Position();
}

//////////////////////////////////////////////////////////////////////

Position ClumpyTorusSpatialGrid::randomPositionInCell(int m) const
{
    if (m >= 0 && m < _numClumps) return randomPositionInSphere(random(), _clumps[m].center(), _clumps[m].radius());
    if (m == _numClumps)
    {
        while (true)
        {
            Position candidate = randomPositionInTorus(random(), _openingAngle, _minRadius, _maxRadius);
            // reject the position if it is inside one of the clumps
            if (cellIndex(candidate) == _numClumps) return candidate;
        }
    }
    return Position();
}

//////////////////////////////////////////////////////////////////////

class ClumpyTorusSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const ClumpyTorusSpatialGrid* _grid{nullptr};

public:
    MySegmentGenerator(const ClumpyTorusSpatialGrid* grid) : _grid(grid) {}

    bool next() override
    {
        switch (state())
        {
            case State::Unknown:
            {
                (void)_grid;
            }

            // intentionally falls through
            case State::Inside:
            {
                (void)_grid;
            }

            // intentionally falls through
            case State::Outside:
            {
                (void)_grid;
            }
        }
        return false;
    }
};

//////////////////////////////////////////////////////////////////////

std::unique_ptr<PathSegmentGenerator> ClumpyTorusSpatialGrid::createPathSegmentGenerator() const
{
    return std::make_unique<MySegmentGenerator>(this);
}

//////////////////////////////////////////////////////////////////////

void ClumpyTorusSpatialGrid::write_xy(SpatialGridPlotFile* outfile) const
{
    (void)outfile;
}

//////////////////////////////////////////////////////////////////////

void ClumpyTorusSpatialGrid::write_xz(SpatialGridPlotFile* outfile) const
{
    (void)outfile;
}

//////////////////////////////////////////////////////////////////////

void ClumpyTorusSpatialGrid::write_yz(SpatialGridPlotFile* outfile) const
{
    (void)outfile;
}

//////////////////////////////////////////////////////////////////////

void ClumpyTorusSpatialGrid::write_xyz(SpatialGridPlotFile* outfile) const
{
    (void)outfile;
}

//////////////////////////////////////////////////////////////////////
