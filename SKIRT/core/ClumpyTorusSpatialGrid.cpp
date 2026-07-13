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
    bool doSpheresOverlap(Position center1, double radius1, Position center2, double radius2)
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

namespace
{
    // Number of candidate split planes ("bins") evaluated per axis during the SAH search
    // (binned SAH, Wald & Havran); O(N) work per node rather than O(N log N) for an exact
    // search, at negligible cost to tree quality.
    constexpr int NumBins = 16;

    // A leaf never holds more clumps than this, regardless of what SAH suggests. Bounds the
    // cost of the linear scan at each leaf, and guarantees the build terminates even when
    // many clumps have coincident or near-coincident centers.
    constexpr int MaxLeafSize = 4;
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

    // Recursively builds the BVH using binned SAH over the clumps referenced by
    // _index[begin,end), reordering that sub-range in place. Appends the newly created
    // node(s) to _nodes and returns the index of the node representing this sub-range.
    // Uses Clump::bounds()/center() directly rather than a cached per-entity box array
    // (unlike the generic BVH class) since neither query here needs individual entity
    // boxes outside of the build itself -- containment and overlap tests below go
    // straight to the exact sphere math instead.
    int buildRecursive(int begin, int end)
    {
        Box box = (*_clumps)[_index[begin]].bounds();
        for (int i = begin + 1; i != end; ++i) box.extend((*_clumps)[_index[i]].bounds());

        int numPrims = end - begin;
        int nodeIndex = static_cast<int>(_nodes.size());

        auto makeLeaf = [&]() {
            _nodes.emplace_back(box, -1, -1, begin, numPrims);
            return nodeIndex;
        };

        if (numPrims <= MaxLeafSize) return makeLeaf();

        // centroid bounding box, used to pick a split axis and the bin boundaries along it
        Position c0 = (*_clumps)[_index[begin]].center();
        double cmin[3] = {c0.x(), c0.y(), c0.z()};
        double cmax[3] = {c0.x(), c0.y(), c0.z()};
        for (int i = begin + 1; i != end; ++i)
        {
            Position c = (*_clumps)[_index[i]].center();
            double v[3] = {c.x(), c.y(), c.z()};
            for (int a = 0; a != 3; ++a)
            {
                cmin[a] = min(cmin[a], v[a]);
                cmax[a] = max(cmax[a], v[a]);
            }
        }
        int axis = 0;
        double extentOnAxis = cmax[0] - cmin[0];
        for (int a = 1; a != 3; ++a)
        {
            double e = cmax[a] - cmin[a];
            if (e > extentOnAxis)
            {
                extentOnAxis = e;
                axis = a;
            }
        }

        // all centroids coincide -> no split can separate the clumps; make a leaf
        // (also protects the binning step below from a division by zero)
        if (extentOnAxis <= 0.0) return makeLeaf();

        auto centroidAxis = [axis](Position c) { return axis == 0 ? c.x() : axis == 1 ? c.y() : c.z(); };

        // bin the clumps of [begin,end) along the chosen axis
        struct Bin
        {
            int count = 0;
            Box box;
            bool hasBox = false;
        };
        vector<Bin> bins(NumBins);
        double lo = cmin[axis], hi = cmax[axis];
        auto binIndexOf = [&](double v) {
            int b = static_cast<int>(NumBins * (v - lo) / (hi - lo));
            return min(max(b, 0), NumBins - 1);
        };
        for (int i = begin; i != end; ++i)
        {
            const Clump& c = (*_clumps)[_index[i]];
            int b = binIndexOf(centroidAxis(c.center()));
            Box eb = c.bounds();
            if (bins[b].hasBox)
                bins[b].box.extend(eb);
            else
            {
                bins[b].box = eb;
                bins[b].hasBox = true;
            }
            bins[b].count++;
        }

        // running unions/counts from the left and from the right, used to evaluate the SAH
        // cost of splitting after each bin boundary without re-scanning for every candidate
        vector<Box> prefixBox(NumBins);
        vector<int> prefixCount(NumBins);
        vector<Box> suffixBox(NumBins);
        vector<int> suffixCount(NumBins);
        {
            bool has = false;
            Box running;
            int runningCount = 0;
            for (int b = 0; b != NumBins; ++b)
            {
                if (bins[b].hasBox)
                {
                    if (has)
                        running.extend(bins[b].box);
                    else
                    {
                        running = bins[b].box;
                        has = true;
                    }
                }
                runningCount += bins[b].count;
                prefixBox[b] = running;
                prefixCount[b] = runningCount;
            }
        }
        {
            bool has = false;
            Box running;
            int runningCount = 0;
            for (int b = NumBins - 1; b >= 0; --b)
            {
                if (bins[b].hasBox)
                {
                    if (has)
                        running.extend(bins[b].box);
                    else
                    {
                        running = bins[b].box;
                        has = true;
                    }
                }
                runningCount += bins[b].count;
                suffixBox[b] = running;
                suffixCount[b] = runningCount;
            }
        }

        // evaluate the SAH cost of splitting right after bin b, keep the cheapest
        double parentArea = box.surfaceArea();
        double bestCost = std::numeric_limits<double>::infinity();
        int bestSplit = -1;
        for (int b = 0; b != NumBins - 1; ++b)
        {
            if (prefixCount[b] == 0 || suffixCount[b + 1] == 0) continue;
            double cost =
                (prefixCount[b] * prefixBox[b].surfaceArea() + suffixCount[b + 1] * suffixBox[b + 1].surfaceArea())
                / parentArea;
            if (cost < bestCost)
            {
                bestCost = cost;
                bestSplit = b;
            }
        }

        // no viable split, or splitting is not cheaper than a leaf -> make a leaf
        double leafCost = static_cast<double>(numPrims);
        if (bestSplit < 0 || bestCost >= leafCost) return makeLeaf();

        double splitPos = lo + (hi - lo) * (bestSplit + 1) / NumBins;
        auto midIt = std::partition(_index.begin() + begin, _index.begin() + end,
                                    [&](int m) { return centroidAxis((*_clumps)[m].center()) < splitPos; });
        int mid = static_cast<int>(midIt - _index.begin());
        if (mid == begin || mid == end) mid = (begin + end) / 2;  // guard against a degenerate partition

        // reserve this node's slot now; children are appended afterwards and get higher
        // indices, so we come back and fill in left/right once the recursion returns
        _nodes.emplace_back(box, -1, -1, 0, 0);
        int leftIndex = buildRecursive(begin, mid);
        int rightIndex = buildRecursive(mid, end);
        _nodes[nodeIndex].left = leftIndex;
        _nodes[nodeIndex].right = rightIndex;
        return nodeIndex;
    }

public:
    BVH() {}

    void loadClumps(const vector<Clump>& clumps)
    {
        _clumps = &clumps;
        int numClumps = static_cast<int>(clumps.size());

        _nodes.clear();
        _index.resize(numClumps);
        for (int m = 0; m != numClumps; ++m) _index[m] = m;

        if (numClumps > 0)
        {
            _nodes.reserve(2 * numClumps / MaxLeafSize + 1);  // rough heuristic upper bound
            buildRecursive(0, numClumps);
        }
    }

    // Returns the indices of a maximal subset of mutually non-overlapping clumps, built
    // greedily in order of increasing index: clump 0 is always kept, and each subsequent
    // clump is kept unless it overlaps a clump already kept. Uses the BVH already built
    // over ALL clumps (including ones that will be rejected here) purely to prune the
    // overlap search down to nearby candidates via box intersection, rather than checking
    // against every previously-accepted clump directly (which would be O(M^2) overall).
    vector<int> allDisjointClumps() const
    {
        vector<int> result;
        int numClumps = static_cast<int>(_index.size());
        if (numClumps == 0) return result;

        vector<char> accepted(numClumps, 0);
        accepted[0] = 1;
        result.push_back(0);

        vector<int> stack;
        for (int i = 1; i != numClumps; ++i)
        {
            const Clump& ci = (*_clumps)[i];
            Box queryBox = ci.bounds();
            bool conflict = false;

            stack.clear();
            stack.push_back(0);
            while (!stack.empty() && !conflict)
            {
                int nodeIndex = stack.back();
                stack.pop_back();
                const Node& node = _nodes[nodeIndex];
                if (!node.box.intersects(queryBox)) continue;

                if (node.isLeaf())
                {
                    for (int k = node.firstIndex; k != node.firstIndex + node.numIndices; ++k)
                    {
                        int m = _index[k];
                        if (!accepted[m]) continue;
                        const Clump& cm = (*_clumps)[m];
                        if (doSpheresOverlap(ci.center(), ci.radius(), cm.center(), cm.radius()))
                        {
                            conflict = true;
                            break;
                        }
                    }
                }
                else
                {
                    stack.push_back(node.left);
                    stack.push_back(node.right);
                }
            }

            if (!conflict)
            {
                accepted[i] = 1;
                result.push_back(i);
            }
        }
        return result;
    }

    // Returns the index of any clump containing the given position, or -1 if none does.
    // No ordering is needed (clumps are non-overlapping by the time this is used in
    // practice), so a plain box-pruned depth-first search returning on the first hit is
    // the fastest option -- no heap, no per-entity box cache, just the exact sphere test
    // at each candidate leaf.
    int anyClumpContaining(Vec bfr) const
    {
        if (_nodes.empty()) return -1;

        thread_local vector<int> stack;
        stack.clear();
        stack.push_back(0);

        while (!stack.empty())
        {
            int nodeIndex = stack.back();
            stack.pop_back();
            const Node& node = _nodes[nodeIndex];
            if (!node.box.contains(bfr)) continue;

            if (node.isLeaf())
            {
                for (int k = node.firstIndex; k != node.firstIndex + node.numIndices; ++k)
                {
                    int m = _index[k];
                    const Clump& c = (*_clumps)[m];
                    if ((bfr - c.center()).norm2() <= square(c.radius())) return m;
                }
            }
            else
            {
                stack.push_back(node.left);
                stack.push_back(node.right);
            }
        }
        return -1;
    }
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
