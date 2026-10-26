/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SelfHeatingRebalancer.hpp"
#include "IterationHistory.hpp"
#include "Log.hpp"
#include "ProcessManager.hpp"
#include "StringUtils.hpp"
#include <array>

////////////////////////////////////////////////////////////////////

namespace
{
    // the ids of the iteration history series recorded by the rebalancer
    enum HistoryId {
        NumRegions,
        NumFallbacks,
        SolveAccepted,
        GlobalFactor,
        MinFactor,
        MaxFactor,
        EffectiveFactor,
        TalliedFraction,
        MinClosure,
        MaxClosure,
        NumSeries
    };

    // the descriptions of the iteration history series recorded by the rebalancer, indexed on their ids
    constexpr std::array<const char*, NumSeries> seriesDescriptions = {
        "number of regions",      "number of regions using the global factor",
        "solve accepted",         "global factor",
        "smallest region factor", "largest region factor",
        "effective factor",       "tallied self-absorbed fraction",
        "smallest tally closure", "largest tally closure"};

    // the range to which the region factors are clamped
    constexpr double minFactor = 0.5;
    constexpr double maxFactor = 5.;

    // the ratio between a pivot and the largest element in its row below which the solution is rejected
    constexpr double minPivotRatio = 1e-12;

    // returns the bin index for each of the specified finite coordinates with nonnegative finite weights, for the
    // specified number of bins with borders at the exact weighted quantiles of the coordinates, so that each bin holds
    // about the same fraction of the total weight; the coordinates are sorted, and each group of equal coordinates
    // goes to the bin holding the cumulative weight at the center of the group, so that equal coordinates always
    // share a bin; if the total weight is zero, all coordinates go to the first bin
    vector<int> quantileBins(const vector<double>& xv, const vector<double>& wv, int numBins)
    {
        size_t n = xv.size();
        vector<int> binv(n, 0);
        if (!n || numBins < 2) return binv;

        double total = 0.;
        for (double w : wv) total += w;
        if (!(total > 0.) || !std::isfinite(total)) return binv;

        // sort the indices of the coordinates in a total order, with ties broken by index, so that the cumulative
        // weights do not depend on the sort implementation
        vector<int> orderv(n);
        std::iota(orderv.begin(), orderv.end(), 0);
        std::sort(orderv.begin(), orderv.end(),
                  [&xv](int a, int b) { return xv[a] < xv[b] || (xv[a] == xv[b] && a < b); });

        // assign each group of equal coordinates to the bin holding the cumulative weight at the center of the group
        double cumulative = 0.;
        for (size_t first = 0; first != n;)
        {
            double x = xv[orderv[first]];
            double weight = wv[orderv[first]];
            size_t last = first + 1;
            for (; last != n && xv[orderv[last]] == x; ++last) weight += wv[orderv[last]];
            int bin = min(numBins - 1, static_cast<int>((cumulative + 0.5 * weight) / total * numBins));
            for (size_t k = first; k != last; ++k) binv[orderv[k]] = bin;
            cumulative += weight;
            first = last;
        }
        return binv;
    }

    // solves the dense linear system a x = b of order n by Gaussian elimination with partial pivoting, where the
    // matrix a is stored in row-major order; on return, b holds the solution and a has been overwritten; returns false
    // if a pivot is smaller than minPivotRatio times the largest element in its row, in which case b is undefined
    bool solveDenseSystem(vector<double>& a, vector<double>& b, size_t n)
    {
        for (size_t k = 0; k != n; ++k)
        {
            // select the pivot row
            size_t p = k;
            for (size_t i = k + 1; i != n; ++i)
                if (abs(a[i * n + k]) > abs(a[p * n + k])) p = i;

            // reject a pivot that is negligible compared to the largest element in its row (this also catches NaNs)
            double rowMax = 0.;
            for (size_t j = k; j != n; ++j) rowMax = max(rowMax, abs(a[p * n + j]));
            if (!(rowMax > 0.) || !(abs(a[p * n + k]) >= minPivotRatio * rowMax)) return false;

            // move the pivot row into place
            if (p != k)
            {
                for (size_t j = 0; j != n; ++j) std::swap(a[k * n + j], a[p * n + j]);
                std::swap(b[k], b[p]);
            }

            // eliminate the column below the pivot
            for (size_t i = k + 1; i != n; ++i)
            {
                double f = a[i * n + k] / a[k * n + k];
                if (f == 0.) continue;
                for (size_t j = k; j != n; ++j) a[i * n + j] -= f * a[k * n + j];
                b[i] -= f * b[k];
            }
        }

        // perform back substitution
        for (size_t i = n; i-- > 0;)
        {
            double sum = b[i];
            for (size_t j = i + 1; j != n; ++j) sum -= a[i * n + j] * b[j];
            b[i] = sum / a[i * n + i];
        }
        return true;
    }

    // returns the specified value clamped to the range of the region factors
    double clampFactor(double f)
    {
        return max(minFactor, min(maxFactor, f));
    }
}

////////////////////////////////////////////////////////////////////

void SelfHeatingRebalancer::Tallies::reset(size_t K)
{
    numRegions = K;
    launched.resize(K);
    absorbed.resize(K * K);
    escaped.resize(K);
}

////////////////////////////////////////////////////////////////////

SelfHeatingRebalancer::SelfHeatingRebalancer(int numAbsorptionRatioBins, int numHeatingRateBins,
                                             IterationHistory* history, const SimulationItem* key, Log* log)
    : _numRatioBins{numAbsorptionRatioBins}, _numHeatingBins{numHeatingRateBins}, _history{history}, _key{key},
      _log{log}
{}

////////////////////////////////////////////////////////////////////

void SelfHeatingRebalancer::begin(const Array& absorbedPrimaryv, const Array& absorbedSecondaryv,
                                  const Array& dustMassv)
{
    int numCells = static_cast<int>(absorbedPrimaryv.size());

    // declare the history series, so that they exist from the first iteration of a loop
    for (int id = 0; id != NumSeries; ++id)
        _history->scalarSeries(_key, id, 1, IterationHistory::Lifetime::Loop, seriesDescriptions[id]);

    // without a secondary radiation field, as in the first iteration of a loop, the regions cannot be formed
    _rebalancing = absorbedSecondaryv.sum() > 0.;
    _numRegions = 0;
    _factorv.clear();
    if (!_rebalancing)
    {
        _active = false;
        return;
    }

    // select the cells that emit and contain dust, and determine the range of the finite ratio coordinates; the
    // luminosity and mass of a selected cell are finite, so that its coordinates are finite as well
    vector<int> cellv;
    double umin = std::numeric_limits<double>::infinity();
    double umax = -umin;
    for (int m = 0; m != numCells; ++m)
    {
        double A = absorbedPrimaryv[m];
        double S = absorbedSecondaryv[m];
        double E = A + S;
        double M = dustMassv[m];
        if (E > 0. && std::isfinite(E) && M > 0. && std::isfinite(M))
        {
            cellv.push_back(m);
            if (A > 0. && S > 0.)
            {
                double u = log10(S) - log10(A);  // a difference of logarithms cannot overflow
                umin = min(umin, u);
                umax = max(umax, u);
            }
        }
    }
    if (!(umin <= umax)) umin = umax = 0.;

    // calculate the coordinates and weights of the selected cells; a cell without absorbed primary luminosity
    // receives the largest finite ratio coordinate, and a cell without absorbed secondary luminosity the smallest one
    size_t n = cellv.size();
    vector<double> uv(n), vv(n), wv(n);
    for (size_t i = 0; i != n; ++i)
    {
        int m = cellv[i];
        double A = absorbedPrimaryv[m];
        double S = absorbedSecondaryv[m];
        double E = A + S;
        uv[i] = A > 0. ? (S > 0. ? log10(S) - log10(A) : umin) : umax;
        vv[i] = log10(E) - log10(dustMassv[m]);
        wv[i] = E;
    }

    // bin each axis separately at the weighted quantiles of its coordinate
    vector<int> iuv = quantileBins(uv, wv, _numRatioBins);
    vector<int> ivv = quantileBins(vv, wv, _numHeatingBins);

    // number the combinations of bins that hold at least one cell consecutively
    vector<int> indexv(_numRatioBins * _numHeatingBins, -1);
    for (size_t i = 0; i != n; ++i) indexv[iuv[i] * _numHeatingBins + ivv[i]] = 0;
    int numRegions = 0;
    for (int& index : indexv)
        if (index == 0) index = numRegions++;
    _numRegions = numRegions;

    // assign each selected cell to its region
    _regionv.assign(numCells, -1);
    for (size_t i = 0; i != n; ++i) _regionv[cellv[i]] = indexv[iuv[i] * _numHeatingBins + ivv[i]];

    // clear the tallies of all execution threads that used them before, and activate the tallies
    for (Tallies* tallies : _tallies.all()) tallies->reset(_numRegions);
    _active = _numRegions > 0;
}

////////////////////////////////////////////////////////////////////

SelfHeatingRebalancer::Tally SelfHeatingRebalancer::tally(int emissionCell)
{
    Tally handle;
    if (_active && emissionCell >= 0)
    {
        int s = _regionv[emissionCell];
        if (s >= 0)
        {
            // a thread that uses the tallies for the first time since begin() receives a new, empty tally object
            Tallies* tallies = _tallies.local();
            if (tallies->numRegions != _numRegions) tallies->reset(_numRegions);

            handle._tallies = tallies;
            handle._regionv = &_regionv;
            handle._numRegions = _numRegions;
            handle._s = s;
        }
    }
    return handle;
}

////////////////////////////////////////////////////////////////////

void SelfHeatingRebalancer::end(const Array& absorbedPrimaryv, const Array& absorbedSecondaryv)
{
    _active = false;
    if (!_rebalancing) return;
    _rebalancing = false;
    size_t K = _numRegions;
    int numCells = static_cast<int>(absorbedPrimaryv.size());

    // records the specified value in the iteration history series with the specified id
    auto record = [this](HistoryId id, double value) {
        _history->scalarSeries(_key, id, 1, IterationHistory::Lifetime::Loop, seriesDescriptions[id]).set(value);
    };

    // without regions, there is nothing to rebalance
    if (!K)
    {
        record(NumRegions, 0.);
        _log->info("  Rebalance: no regions, because no cell with dust emits");
        return;
    }

    // sum the tallies over execution threads and processes, combined in a single array to limit communication:
    // the absorbed luminosities at indices r*K+s, followed by the launched and escaped luminosities
    Array sumv(K * K + 2 * K);
    for (const Tallies* tallies : _tallies.all())
    {
        if (tallies->numRegions != K) continue;
        for (size_t i = 0; i != K * K; ++i) sumv[i] += tallies->absorbed[i];
        for (size_t s = 0; s != K; ++s)
        {
            sumv[K * K + s] += tallies->launched[s];
            sumv[K * K + K + s] += tallies->escaped[s];
        }
    }
    ProcessManager::sumToAll(sumv);
    auto absorbed = [&sumv, K](size_t r, size_t s) { return sumv[r * K + s]; };
    auto launched = [&sumv, K](size_t s) { return sumv[K * K + s]; };
    auto escaped = [&sumv, K](size_t s) { return sumv[K * K + K + s]; };

    // calculate the absorbed primary and secondary luminosity of each region, and the totals over all cells
    vector<double> Av(K, 0.), Cv(K, 0.);
    double totalSecondary = 0.;
    for (int m = 0; m != numCells; ++m)
    {
        int r = _regionv[m];
        if (r >= 0)
        {
            Av[r] += absorbedPrimaryv[m];
            Cv[r] += absorbedSecondaryv[m];
        }
        totalSecondary += absorbedSecondaryv[m];
    }
    double sumA = 0., sumC = 0., sumLaunched = 0., sumAbsorbed = 0.;
    for (size_t r = 0; r != K; ++r)
    {
        sumA += Av[r];
        sumC += Cv[r];
        sumLaunched += launched(r);
        for (size_t s = 0; s != K; ++s) sumAbsorbed += absorbed(r, s);
    }

    // calculate the global factor, i.e. the factor of a rebalance with a single region
    double p = sumLaunched > 0. ? sumAbsorbed / sumLaunched : 0.;
    double globalFactor = (sumLaunched > 0. && p < 1. && sumC > 0.) ? (sumA / (1. - p) - sumA) / sumC : 1.;
    globalFactor = clampFactor(globalFactor);

    // set up the system (I - T) x = A, where T holds the transfer fractions, and solve it
    vector<double> matrixv(K * K, 0.);
    for (size_t s = 0; s != K; ++s)
    {
        double Ls = launched(s);
        for (size_t r = 0; r != K; ++r) matrixv[r * K + s] = (r == s ? 1. : 0.) - (Ls > 0. ? absorbed(r, s) / Ls : 0.);
    }
    vector<double> xv(Av);
    bool accepted = solveDenseSystem(matrixv, xv, K);
    if (accepted)
        for (double x : xv)
            if (!(x >= 0.)) accepted = false;

    // determine the factor of each region, falling back to the global factor where needed
    _factorv.assign(K, globalFactor);
    int numFallbacks = 0;
    for (size_t r = 0; r != K; ++r)
    {
        if (accepted && launched(r) > 0. && Cv[r] > 0. && xv[r] >= Av[r])
            _factorv[r] = clampFactor((xv[r] - Av[r]) / Cv[r]);
        else
            numFallbacks++;
    }
    auto [minIt, maxIt] = std::minmax_element(_factorv.begin(), _factorv.end());
    double minRegionFactor = *minIt;
    double maxRegionFactor = *maxIt;

    // calculate the effective factor, i.e. the ratio of the total absorbed secondary luminosity after and before the
    // rebalance, taking into account that cells without a region keep their value
    double rescaledSecondary = totalSecondary - sumC;
    for (size_t r = 0; r != K; ++r) rescaledSecondary += _factorv[r] * Cv[r];
    double effectiveFactor = totalSecondary > 0. ? rescaledSecondary / totalSecondary : 1.;

    // determine the range of the tally closure over the regions from which luminosity was launched
    double minClosure = std::numeric_limits<double>::infinity();
    double maxClosure = -minClosure;
    for (size_t s = 0; s != K; ++s)
    {
        if (launched(s) > 0.)
        {
            double sum = escaped(s);
            for (size_t r = 0; r != K; ++r) sum += absorbed(r, s);
            minClosure = min(minClosure, sum / launched(s));
            maxClosure = max(maxClosure, sum / launched(s));
        }
    }

    // record the rebalance series in the iteration history
    record(NumRegions, static_cast<double>(K));
    record(NumFallbacks, numFallbacks);
    record(SolveAccepted, accepted ? 1. : 0.);
    record(GlobalFactor, globalFactor);
    record(MinFactor, minRegionFactor);
    record(MaxFactor, maxRegionFactor);
    record(EffectiveFactor, effectiveFactor);
    record(TalliedFraction, p);
    if (minClosure <= maxClosure)
    {
        record(MinClosure, minClosure);
        record(MaxClosure, maxClosure);
    }

    // log a summary
    string message = "  Rebalance: " + std::to_string(K) + " regions, ";
    if (accepted)
    {
        message += "factors " + StringUtils::toString(minRegionFactor, 'f', 3) + " to "
                   + StringUtils::toString(maxRegionFactor, 'f', 3) + ", effective factor "
                   + StringUtils::toString(effectiveFactor, 'f', 4) + ", ";
        if (numFallbacks)
            message += std::to_string(numFallbacks) + " using the global factor "
                       + StringUtils::toString(globalFactor, 'f', 4);
        else
            message += "no fallbacks";
    }
    else
    {
        message +=
            "solve rejected, global factor " + StringUtils::toString(globalFactor, 'f', 4) + " applied to all regions";
    }
    _log->info(message);
}

////////////////////////////////////////////////////////////////////

double SelfHeatingRebalancer::factor(int m) const
{
    int r = _regionv.empty() ? -1 : _regionv[m];
    return r >= 0 && !_factorv.empty() ? _factorv[r] : 1.;
}

////////////////////////////////////////////////////////////////////
