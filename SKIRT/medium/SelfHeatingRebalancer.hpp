/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SELFHEATINGREBALANCER_HPP
#define SELFHEATINGREBALANCER_HPP

#include "Array.hpp"
#include "ThreadLocalMember.hpp"
class IterationHistory;
class Log;
class SimulationItem;

////////////////////////////////////////////////////////////////////

/** SelfHeatingRebalancer accelerates the iteration over secondary emission in a simulation in which
    thermal dust emission is the only form of secondary emission, by rebalancing the secondary
    radiation field over regions of spatial cells after each iteration. The medium system creates
    and owns a rebalancer if the simulation iterates over secondary emission and the user did not
    disable the rebalance.

    Without acceleration, the iteration over secondary emission is a Lambda iteration: its error
    shrinks per iteration by a factor close to the fraction of the dust emission that is reabsorbed
    by the dust, so that strongly self-absorbing models need very many iterations. The rebalance
    adapts the coarse-mesh rebalance used for Monte Carlo neutron transport. During each secondary
    emission iteration, it records how much of the dust emission launched from each region is
    absorbed by the dust in each other region, and how much escapes. At the end of the iteration, it
    solves a small linear system for the dust emission of each region that is consistent with these
    transfer fractions and with the absorbed primary luminosity of each region, and it determines a
    factor for each region by which the medium system rescales the stable secondary radiation field
    and the stable absorbed secondary luminosity of all cells in the region. The next iteration
    then starts from that luminosity, while the spectral shape of the radiation field in each cell
    is kept. At the fixed point, all factors equal one, so that the rebalance does not change the
    converged solution.

    <b>Regions</b>

    A region is not a contiguous part of the spatial domain: it collects cells in similar physical
    conditions, which may be scattered throughout the domain. The regions are formed from the
    absorbed luminosities that determine the emission of the current iteration: the absorbed
    primary luminosity \f$A_m\f$ and the stable absorbed secondary luminosity \f$S_m\f$ of each
    cell \f$m\f$, which add up to the emitted luminosity \f$E_m=A_m+S_m\f$, and the dust mass
    \f$M_m\f$ of the cell. Cells with \f$E_m=0\f$ or \f$M_m=0\f$ belong to no region. For the other
    cells, the rebalancer calculates the coordinates \f$u_m=\log_{10}(S_m/A_m)\f$, which measures
    how deeply a cell sits in the self-absorbing core, and \f$v_m=\log_{10}(E_m/M_m)\f$, which
    measures the heating. A cell without absorbed primary luminosity receives the largest finite
    value of \f$u\f$ over all cells, and a cell without absorbed secondary luminosity the smallest
    one. Each axis is divided into bins with borders at the weighted quantiles of the coordinate,
    with the emitted luminosity as weight, so that each bin holds about the same fraction of the
    emission. The quantiles are exact: the cells are sorted on the coordinate, and each group of
    cells with the same coordinate goes to the bin holding the cumulative emission at the center of
    the group, so that such cells, for example all cells without absorbed secondary luminosity,
    always share a bin. Unlike a histogram, this does not depend on the range of the coordinates,
    which a few cells with a tiny absorbed luminosity can stretch over hundreds of decades. The
    regions are the combinations of a bin on each axis that hold at least one cell, numbered
    consecutively.

    <b>Tallies</b>

    During the secondary emission segment of an iteration, the medium system obtains a Tally handle
    for each path of a dust-emitted photon packet, and adds through it the launched luminosity (for
    the first path of the packet), the luminosity absorbed by the dust in each cell crossed by the
    path, and the luminosity leaving the spatial domain at the end of the path. The handle binds the
    tally object of the current execution thread to the region of the cell from which the packet
    was emitted. Each thread accumulates into its own tally object, managed through a
    ThreadLocalMember, to avoid contention. The begin() and end() functions, which are called
    between segments, zero and sum these objects, respectively.

    <b>Rebalance system</b>

    With the tallies summed over threads and processes, the end() function calculates the transfer
    fractions \f$T_{rs}\f$, i.e. the fraction of the luminosity launched from region \f$s\f$ that is
    absorbed in region \f$r\f$, the absorbed primary luminosity \f$A_r\f$ of each region, and the
    absorbed secondary luminosity \f$C_r\f$ of each region in the current iteration. It then solves
    the system \f$(I-T)\,x=A\f$ for the emission \f$x_r\f$ of each region by Gaussian elimination
    with partial pivoting. The factor of region \f$r\f$ is \f$f_r=(x_r-A_r)/C_r\f$.

    The following safeguards protect against noisy or degenerate tallies. The global factor is the
    factor of a rebalance with a single region, calculated from the tallied self-absorbed fraction
    \f$p\f$ (the absorbed luminosity divided by the launched luminosity, over all regions) as
    \f$\left(\sum A_r/(1-p) - \sum A_r\right) / \sum C_r\f$, or one if \f$p\ge1\f$ or \f$\sum
    C_r=0\f$. A region receives the global factor if no luminosity was launched from it, if it has
    no absorbed secondary luminosity, or if \f$x_r<A_r\f$. If a pivot is smaller than
    \f$10^{-12}\f$ times the largest element of its row, or if any \f$x_r\f$ is negative, the
    solution is rejected and all regions receive the global factor. Finally, all factors are
    clamped to the range from 0.5 to 5.

    The end() function logs a summary of the rebalance and records the corresponding series in the
    iteration history. These series include the smallest and largest tally closure over the
    regions, i.e. the absorbed plus escaped luminosity divided by the launched luminosity. Apart
    from absorption by other media and the luminosity lost when photon packets are terminated, the
    closure equals one within Monte Carlo noise.

    The partition, the solution, and the factors are calculated identically on all processes,
    because their inputs are identical on all processes. */
class SelfHeatingRebalancer
{
    //============= Tallies =============

private:
    // the tallies accumulated by a single execution thread
    struct Tallies
    {
        size_t numRegions{0};  // the number of regions for which the tallies are sized
        Array launched;        // the luminosity launched from each region s
        Array absorbed;        // the luminosity launched from region s and absorbed in region r, indexed on r*K+s
        Array escaped;         // the luminosity launched from each region s that leaves the spatial domain

        // resizes the tallies for the specified number of regions and sets all values to zero
        void reset(size_t K);
    };

public:
    /** A Tally instance is a handle through which the medium system adds the contributions of a
        single photon packet path to the tallies of the current execution thread. A default
        constructed instance is inactive and converts to false. An active instance is obtained from
        the tally() function and converts to true. The functions adding contributions may be called
        only for an active instance. A Tally instance is valid only during the secondary emission
        segment in which it was obtained. */
    class Tally
    {
    public:
        /** This function returns true if the tally is active, and false otherwise. */
        explicit operator bool() const { return _tallies != nullptr; }

        /** This function adds the specified luminosity to the luminosity launched from the region
            of the emission cell. It should be called for the first path of a photon packet only. */
        void launch(double L) const { _tallies->launched[_s] += L; }

        /** This function adds the specified luminosity to the luminosity launched from the region
            of the emission cell and absorbed in the region of the spatial cell with index \em m.
            If that cell belongs to no region, the function does nothing. */
        void absorb(int m, double L) const
        {
            int r = (*_regionv)[m];
            if (r >= 0) _tallies->absorbed[static_cast<size_t>(r) * _numRegions + _s] += L;
        }

        /** This function adds the specified luminosity to the luminosity launched from the region
            of the emission cell that leaves the spatial domain. */
        void escape(double L) const { _tallies->escaped[_s] += L; }

    private:
        friend class SelfHeatingRebalancer;
        Tallies* _tallies{nullptr};
        const vector<int>* _regionv{nullptr};
        size_t _numRegions{0};
        size_t _s{0};
    };

    //============= Construction =============

public:
    /** The constructor creates a rebalancer with the specified number of bins on the axis of the
        ratio of the absorbed secondary to the absorbed primary luminosity and on the axis of the
        absorbed luminosity per unit dust mass. It also receives the iteration history in which the
        rebalance series are recorded, the simulation item on which these series are keyed, and the
        logger for the summary messages. */
    SelfHeatingRebalancer(int numAbsorptionRatioBins, int numHeatingRateBins, IterationHistory* history,
                          const SimulationItem* key, Log* log);

    //============= Rebalancing =============

public:
    /** This function partitions the spatial cells into regions as described in the class header,
        based on the specified absorbed primary luminosity, stable absorbed secondary luminosity,
        and dust mass of each cell, clears the tallies of all execution threads, and activates the
        tallies. It must be called before the secondary photon packets of each iteration of a
        secondary or merged iteration loop are launched.

        If the stable absorbed secondary luminosity is zero in all cells, as is the case in the
        first iteration of a loop, the regions cannot be formed. The function then leaves the
        tallies inactive, and the subsequent call to end() does nothing. In any case, the function
        declares the rebalance series in the iteration history, so that they exist from the first
        iteration of the loop, as required by the HistoryProbe. */
    void begin(const vector<double>& absorbedPrimaryv, const vector<double>& absorbedSecondaryv,
               const vector<double>& dustMassv);

    /** This function returns a tally handle for a photon packet path, bound to the tally object of
        the current execution thread and to the region of the specified emission cell. The handle is
        inactive if the tallies are not active, or if the emission cell is negative or belongs to no
        region. The function may be called from multiple execution threads in parallel. */
    Tally tally(int emissionCell);

    /** Unless the corresponding call to begin() found no secondary radiation field, in which case
        it does nothing, this function deactivates the tallies, sums them over execution threads
        and processes, solves the rebalance system, and determines the factor of each region as
        described in the class header, based on the specified absorbed primary luminosity and
        absorbed secondary luminosity of each cell in the current iteration. It logs a summary and
        records the rebalance series in the iteration history. After this function returns, the
        factor() function returns the factor for each cell, which is one if nothing was
        rebalanced. The function must be called after the radiation field of the secondary emission
        segment has been communicated between processes, and before the stable secondary radiation
        field is rescaled. */
    void end(const vector<double>& absorbedPrimaryv, const vector<double>& absorbedSecondaryv);

    /** This function returns the factor determined by the most recent call to end() for the
        spatial cell with index \em m, or one if the cell belongs to no region or if nothing was
        rebalanced. */
    double factor(int m) const;

    //============= Data members =============

private:
    // configuration
    int _numRatioBins{0};
    int _numHeatingBins{0};
    IterationHistory* _history{nullptr};
    const SimulationItem* _key{nullptr};
    Log* _log{nullptr};

    // partition and factors
    bool _rebalancing{false};  // true between begin() and end() if there is a secondary radiation field
    bool _active{false};       // true between begin() and end() if, in addition, there is at least one region
    size_t _numRegions{0};     // the number of regions K
    vector<int> _regionv;      // the region index for each spatial cell, or -1
    vector<double> _factorv;   // the factor for each region, determined by end()

    // tallies for each execution thread
    ThreadLocalMember<Tallies> _tallies;
};

////////////////////////////////////////////////////////////////////

#endif
