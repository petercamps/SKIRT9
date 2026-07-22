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
#include <array>

////////////////////////////////////////////////////////////////////

// Geometric helper functions
namespace
{
    // small value
    constexpr double EPS = 1e-12;

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

    // This function returns true if the given position is inside the given sphere
    bool isPositionInSphere(Vec p, Vec center, double r)
    {
        return (p - center).norm2() <= square(r);
    }

    // This function returns true if the specified spheres overlap or touch, and false otherwise
    bool doSpheresOverlap(Vec center1, double radius1, Vec center2, double radius2)
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

    // Solves a*t^2 + 2*b*t + c = 0 for t, returning the number of distinct real roots found
    // (0, 1, or 2) and storing them in t0 <= t1. Falls back to the linear equation 2*b*t + c = 0
    // when a is (numerically) zero, which happens for the wedge-cone quadratic below whenever the
    // ray direction makes exactly the wedge's half-angle with the equatorial plane.
    int solveQuadratic(double a, double b, double c, double& t0, double& t1)
    {
        if (std::fabs(a) < EPS)
        {
            if (std::fabs(b) < EPS) return 0;
            t0 = -0.5 * c / b;
            return 1;
        }
        double disc = b * b - a * c;
        if (disc < 0.) return 0;
        double sq = std::sqrt(disc);
        double r0 = (-b - sq) / a;
        double r1 = (-b + sq) / a;
        if (r0 <= r1)
        {
            t0 = r0;
            t1 = r1;
        }
        else
        {
            t0 = r1;
            t1 = r0;
        }
        return 2;
    }

    // Computes the intersection(s) of the ray (r0, k) -- k assumed a unit vector -- with the
    // sphere of given center and radius. Returns the number of distinct real roots (0, 1, or 2),
    // stored in t0 <= t1.
    int intersectSphere(Position r0, Direction k, Position center, double radius, double& t0, double& t1)
    {
        Vec d = r0 - center;
        return solveQuadratic(1., Vec::dot(d, k), d.norm2() - square(radius), t0, t1);
    }

    // Computes the intersection(s) of the ray (r0, k) with the double cone that bounds the torus's
    // angular wedge, i.e. the surface |z| = s*tan(delta) where s is the cylindrical radius,
    // expressed in the numerically friendlier factored form
    //     z^2*cos^2(delta) - (x^2+y^2)*sin^2(delta) = 0
    // which avoids dividing by s (a problem on the z-axis) or by cos(delta) (a problem as delta
    // approaches 90 deg). Returns the number of distinct real roots (0, 1, or 2), in t0 <= t1.
    int intersectWedgeCone(Position r0, Direction k, double delta, double& t0, double& t1)
    {
        double c2 = square(std::cos(delta));
        double s2 = square(std::sin(delta));
        double a = c2 * square(k.z()) - s2 * (square(k.x()) + square(k.y()));
        double b = c2 * r0.z() * k.z() - s2 * (r0.x() * k.x() + r0.y() * k.y());
        double cc = c2 * square(r0.z()) - s2 * (square(r0.x()) + square(r0.y()));
        return solveQuadratic(a, b, cc, t0, t1);
    }

    // Returns the smallest strictly positive root of the ray-sphere intersection, or 0 if there is
    // none (either no real intersection, or both intersections lie behind the ray's origin).
    double firstIntersectionSphere(Position r0, Direction k, Position center, double radius)
    {
        double t0, t1;
        int n = intersectSphere(r0, k, center, radius, t0, t1);
        if (n >= 1 && t0 > 0.) return t0;
        if (n == 2 && t1 > 0.) return t1;
        return 0.;
    }

    // Returns the smallest strictly positive root of the ray-wedge-cone intersection, or 0 if
    // there is none.
    double firstIntersectionWedgeCone(Position r0, Direction k, double delta)
    {
        double t0, t1;
        int n = intersectWedgeCone(r0, k, delta, t0, t1);
        if (n >= 1 && t0 > 0.) return t0;
        if (n == 2 && t1 > 0.) return t1;
        return 0.;
    }

    // Given a position known to be INSIDE the torus domain, returns the distance to the nearest
    // domain-boundary crossing ahead along the ray (r0, k). Because all three domain conditions
    // (outside r1, inside r2, inside the wedge) currently hold, crossing ANY one of the three
    // boundary surfaces first is necessarily the moment the domain is exited -- so this is simply
    // the smallest of the three nearest forward crossings. Returns 0 only in the (should-not-occur
    // for a bounded domain) case that no exit is found, handled defensively.
    double distanceToDomainExit(Position r0, Direction k, double delta, double r1, double r2)
    {
        double best = DBL_MAX;
        double t;
        t = firstIntersectionSphere(r0, k, Position(), r1);
        if (t > 0. && t < best) best = t;
        t = firstIntersectionSphere(r0, k, Position(), r2);
        if (t > 0. && t < best) best = t;
        t = firstIntersectionWedgeCone(r0, k, delta);
        if (t > 0. && t < best) best = t;
        return best < DBL_MAX ? best : 0.;
    }

    // Given a position possibly OUTSIDE the torus domain, returns the distance to the first
    // position ahead along the ray (r0, k) where the path (re-)enters the domain, or 0 if it
    // never does.
    //
    // The domain's boundary is only easy to reason about piecewise, since it's built from three
    // separate quadric surfaces. This function collects every forward crossing of any of the
    // three surfaces, sorts them by distance, and sweeps through them while tracking which of the
    // three domain conditions currently holds (each crossing of a surface flips exactly the one
    // condition associated with it). The first crossing after which all three hold simultaneously
    // is the entry point. This handles the torus's non-convexity -- and any number of dips in and
    // out -- correctly and generally, without needing to special-case the geometry.
    double firstDomainEntry(Position r0, Direction k, double delta, double r1, double r2)
    {
        double rho2 = r0.norm2();
        bool outsideInner = rho2 > square(r1);
        bool insideOuter = rho2 < square(r2);
        double sxy2 = square(r0.x()) + square(r0.y());
        bool insideWedge = square(r0.z()) < sxy2 * square(std::tan(delta));

        struct Crossing
        {
            double t;
            int tag;  // 0: crosses r1, 1: crosses r2, 2: crosses the wedge cone
        };

        // at most 2 roots per surface (r1 sphere, r2 sphere, wedge cone) -> 6 crossings max;
        // a fixed-size array avoids a heap allocation in this per-ray hot path. 'count' tracks
        // how many of the 6 slots are actually filled; only crossings[0,count) is ever read.
        std::array<Crossing, 6> crossings;
        int count = 0;

        double t0, t1;
        int n;

        n = intersectSphere(r0, k, Position(), r1, t0, t1);
        if (n >= 1 && t0 > 0.) crossings[count++] = {t0, 0};
        if (n == 2 && t1 > 0.) crossings[count++] = {t1, 0};

        n = intersectSphere(r0, k, Position(), r2, t0, t1);
        if (n >= 1 && t0 > 0.) crossings[count++] = {t0, 1};
        if (n == 2 && t1 > 0.) crossings[count++] = {t1, 1};

        n = intersectWedgeCone(r0, k, delta, t0, t1);
        if (n >= 1 && t0 > 0.) crossings[count++] = {t0, 2};
        if (n == 2 && t1 > 0.) crossings[count++] = {t1, 2};

        std::sort(crossings.begin(), crossings.begin() + count,
                  [](const Crossing& a, const Crossing& b) { return a.t < b.t; });

        for (int i = 0; i != count; ++i)
        {
            const Crossing& x = crossings[i];
            if (x.tag == 0)
                outsideInner = !outsideInner;
            else if (x.tag == 1)
                insideOuter = !insideOuter;
            else
                insideWedge = !insideWedge;

            if (outsideInner && insideOuter && insideWedge) return x.t;
        }
        return 0.;
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

// Private utility class for organizing spheres in a data structure that allows efficient queries.
// The implementation employs a linearized Bounding Volume Hierarchy (BVH) that is bulk-loaded
// using a surface area heuristic (SAH) for optimal balance. Construction of the data structure
// runs in a single serial thread. Because the data structure does not change after initial
// construction, all queries are thread-safe.
class ClumpyTorusSpatialGrid::BVH
{
private:
    class Node
    {
    public:
        Box box;         // union of the bounding boxes of all clumps in this node's subtree
        int left;        // index into _nodes of the left child, or -1 if this is a leaf
        int right;       // index into _nodes of the right child, or -1 if this is a leaf
        int firstIndex;  // for a leaf: start of this leaf's range in _index
        int numIndices;  // for a leaf: number of clumps in this leaf; 0 for an interior node

        Node(const Box& box_, int left_, int right_, int firstIndex_, int numIndices_)
            : box(box_), left(left_), right(right_), firstIndex(firstIndex_), numIndices(numIndices_)
        {}

        bool isLeaf() const { return numIndices > 0; }
    };

    // flat array of BVH nodes; index 0 is the root, meaningful only when _numEntities > 0
    vector<Node> _nodes;

    // entity indices, reordered during the build so that the clumps of any single leaf
    // occupy the contiguous range [firstIndex, firstIndex+numIndices) of this array
    vector<int> _index;

    // pointer to array with clumps, passed on when bulk-loading the BVH
    const vector<Clump>* _clumps;

    // Recursively builds the BVH using binned SAH over the clumps referenced by
    // _index[begin,end), reordering that sub-range in place. Appends the newly created
    // node(s) to _nodes and returns the index of the node representing this sub-range.
    int buildRecursive(int begin, int end)
    {
        // aggregate bounding box for this sub-range: the union of the entity boxes referenced
        // by _index[begin,end).
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

        // evaluate the SAH cost of splitting right after bin b, for every interior bin boundary,
        // and keep the cheapest; cost is expressed relative to the parent's surface area so it
        // can be compared directly against the cost of not splitting at all
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
    // Bulk-loads the specified clumps into the BVH
    void loadClumps(const vector<Clump>& clumps)
    {
        _clumps = &clumps;
        int numClumps = static_cast<int>(clumps.size());

        _nodes.clear();
        _index.resize(numClumps);
        for (int m = 0; m != numClumps; ++m) _index[m] = m;

        if (numClumps > 0)
        {
            // rough heuristic upper bound on the final node count (a binary tree with roughly
            // numEntities/MaxLeafSize leaves has on the order of twice as many nodes total)
            _nodes.reserve(2 * numClumps / MaxLeafSize + 1);
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

        // accepted[m] becomes true once entity m has been added to the result; a plain
        // vector<char> is used instead of vector<bool> for fast, unambiguous random-access.
        vector<char> accepted(numClumps, 0);
        accepted[0] = 1;
        result.push_back(0);

        // For each remaining clump, we need to know whether it overlaps ANY already-accepted clump
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
    // the fastest option.
    int anyClumpContaining(Vec bfr) const
    {
        if (_nodes.empty()) return -1;

        // use a thread_local stack to avoid allocations between consecutibe queries
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
                    if (isPositionInSphere(bfr, c.center(), c.radius())) return m;
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

    // Returns the index of the nearest clump intersected by the ray (r0,k) at a forward
    // distance strictly between 0 and tMax, or -1 if there is none. On success, tEntry and
    // tExit are set to the entry and exit distances of that specific clump.
    //
    // This is a nearest-first, depth-first traversal with branch-and-bound pruning: a node's
    // ray-box entry distance (from Box::intersects) is a valid lower bound on the entry
    // distance of any clump inside it, so a node whose entry distance is not smaller than the
    // best confirmed hit so far can be skipped outright. Visiting the geometrically nearer
    // child first tends to tighten that bound quickly. Unlike an ordered multi-hit iterator,
    // this only ever needs the single nearest hit, since the segment generator re-queries
    // fresh from its new position after every segment.
    int nearestClumpAlongRay(Position r0, Direction k, double tMax, double& tEntry, double& tExit) const
    {
        if (_nodes.empty()) return -1;

        double best = tMax;
        int bestIndex = -1;
        double bestEntry = 0.;
        double bestExit = 0.;

        struct StackEntry
        {
            int nodeIndex;
            double boxEntry;
        };
        thread_local vector<StackEntry> stack;  // thread_local to avoid reallocations on repeated use
        stack.clear();

        double smin, smax;
        if (_nodes[0].box.intersects(r0, k, smin, smax) && smin < best) stack.push_back({0, smin});

        while (!stack.empty())
        {
            StackEntry top = stack.back();
            stack.pop_back();
            if (top.boxEntry >= best) continue;  // superseded by a tighter bound found meanwhile

            const Node& node = _nodes[top.nodeIndex];
            if (node.isLeaf())
            {
                for (int i = node.firstIndex; i != node.firstIndex + node.numIndices; ++i)
                {
                    int m = _index[i];
                    const Clump& c = (*_clumps)[m];
                    double t0, t1;
                    int n = intersectSphere(r0, k, c.center(), c.radius(), t0, t1);
                    if (n >= 1 && t0 > 0. && t0 < best)
                    {
                        best = t0;
                        bestIndex = m;
                        bestEntry = t0;
                        bestExit = (n == 2) ? t1 : t0;
                    }
                }
            }
            else
            {
                double sminL, smaxL, sminR, smaxR;
                bool hitL = _nodes[node.left].box.intersects(r0, k, sminL, smaxL) && sminL < best;
                bool hitR = _nodes[node.right].box.intersects(r0, k, sminR, smaxR) && sminR < best;

                // push the farther child first so the nearer one is popped (and processed) first
                if (hitL && hitR)
                {
                    if (sminL <= sminR)
                    {
                        stack.push_back({node.right, sminR});
                        stack.push_back({node.left, sminL});
                    }
                    else
                    {
                        stack.push_back({node.left, sminL});
                        stack.push_back({node.right, sminR});
                    }
                }
                else if (hitL)
                    stack.push_back({node.left, sminL});
                else if (hitR)
                    stack.push_back({node.right, sminR});
            }
        }

        tEntry = bestEntry;
        tExit = bestExit;
        return bestIndex;
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

    // limit the epsilon we use for progressing the path to a value smaller than the hole and the smallest clump
    double smallestRadius = _minRadius;
    for (const Clump& clump : _clumps)
    {
        double r = clump.radius();
        if (r < smallestRadius) smallestRadius = r;
    }
    _eps = std::min(EPS * _maxRadius, 0.1 * smallestRadius);
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
            // accept the position if it is not inside one of the clumps
            int m = _bvh->anyClumpContaining(candidate);
            if (m < 0) return candidate;
        }
    }
    return Position();
}

//////////////////////////////////////////////////////////////////////

class ClumpyTorusSpatialGrid::MySegmentGenerator : public PathSegmentGenerator
{
    const ClumpyTorusSpatialGrid* _grid{nullptr};

    // the clump index containing the current position, or _grid->_numClumps for the
    // background cell; meaningful only while state() is State::Inside
    int _cell{-1};

public:
    MySegmentGenerator(const ClumpyTorusSpatialGrid* grid) : _grid(grid) {}

    bool next() override
    {
        switch (state())
        {
            case State::Unknown:
            {
                // classify the starting position, then dispatch to whichever of the two
                // cases below actually applies -- fallthrough can't do this conditionally,
                // since it always proceeds to the next case in source order regardless of
                // the state we just set, so we recurse instead
                if (isPositionInTorus(r(), _grid->_openingAngle, _grid->_minRadius, _grid->_maxRadius))
                {
                    _cell = _grid->_bvh->anyClumpContaining(r());
                    if (_cell < 0) _cell = _grid->_numClumps;
                    setState(State::Inside);
                }
                else
                {
                    setState(State::Outside);
                }
                return next();
            }

            case State::Outside:
            {
                // find the next point, if any, where the path (re-)enters the domain; this
                // handles both "never entered yet" and "exited and might come back", since
                // the torus is not convex and a path can dip in and out multiple times
                double ds = firstDomainEntry(r(), k(), _grid->_openingAngle, _grid->_minRadius, _grid->_maxRadius);
                if (ds <= 0.) return false;  // never (re-)enters the domain -- generation is complete

                setEmptySegment(ds);
                propagater(ds + _grid->_eps);

                _cell = _grid->_bvh->anyClumpContaining(r());
                if (_cell < 0) _cell = _grid->_numClumps;
                setState(State::Inside);
                return true;
            }

            case State::Inside:
            {
                double ds;
                if (_cell < _grid->_numClumps)
                {
                    // inside a specific clump: the segment ends where the path exits that
                    // sphere (a single, unambiguous positive root, since the origin is
                    // strictly inside the sphere)
                    const Clump& c = _grid->_clumps[_cell];
                    ds = firstIntersectionSphere(r(), k(), c.center(), c.radius());
                    setSegment(_cell, ds);
                    propagater(ds + _grid->_eps);

                    // clumps are mutually disjoint and don't touch the domain boundary, so
                    // exiting one always lands unambiguously back in the background cell
                    _cell = _grid->_numClumps;
                }
                else
                {
                    // inside the background: the segment ends at whichever comes first, the
                    // domain boundary or the entry point of the nearest upcoming clump; bounding
                    // the BVH search by the domain-exit distance both prunes the search and
                    // means "no clump found" and "found beyond the domain boundary" collapse
                    // into the same outcome
                    double dsDomain =
                        distanceToDomainExit(r(), k(), _grid->_openingAngle, _grid->_minRadius, _grid->_maxRadius);
                    double tEntry, tExit;
                    int m = _grid->_bvh->nearestClumpAlongRay(r(), k(), dsDomain, tEntry, tExit);

                    ds = (m >= 0) ? tEntry : dsDomain;
                    setSegment(_grid->_numClumps, ds);
                    propagater(ds + _grid->_eps);

                    if (m >= 0)
                        _cell = m;
                    else
                        setState(State::Outside);
                }
                return true;
            }
        }
        return false;  // unreachable; silences a compiler warning about a missing return
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
