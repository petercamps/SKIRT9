/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef HISTORYPROBE_HPP
#define HISTORYPROBE_HPP

#include "IterationHistory.hpp"
#include "Probe.hpp"
#include "TextOutFile.hpp"

////////////////////////////////////////////////////////////////////

/** HistoryProbe writes the convergence history of each iteration loop of the simulation to a text
    column file, with a row for each iteration and a column for each scalar series held by the
    IterationHistory. A simulation that iterates over primary emission, over secondary emission, or
    over merged primary and secondary emission thus produces a file for each of these loops, named
    <tt>prefix_probe_primary.dat</tt>, <tt>prefix_probe_secondary.dat</tt>, or
    <tt>prefix_probe_merged.dat</tt>, where \em prefix is the simulation prefix and \em probe is the
    name of the probe. The file gives a compact, machine-readable record of how each loop
    converged, for example to find out which criterion keeps a loop from converging, to tune the
    iteration parameters, or to plot the convergence behavior.

    The first column holds the one-based iteration index within the loop. The other columns follow
    the scalar series that exist in the history when the probe is first performed in the loop, in
    the order in which they were declared. These include, for example, the number of photon packets
    launched in each iteration, whether the loop has converged (0 or 1), the dust-absorbed
    secondary luminosity, and the series kept by material mixes and dynamic state recipes, such as
    aggregates of medium state variables. Each column description combines the description of the
    series with a description of the simulation item on whose behalf it is kept: its type, or for a
    medium component, its index. Values are converted to output units according to the quantity of
    each series. A value that was not set in a given iteration is written as NaN. A series with
    simulation lifetime, such as an aggregate, appears in the file of each loop, so that its
    evolution across loops can be followed by combining the files. A series that is first declared
    after the probe has been performed in a loop is left out of the file for that loop, and the
    probe logs a warning.

    The probe appends a row to the file after each iteration and flushes it, so that the progress
    of a long run can be followed while it executes, and so that the output survives an aborted
    run. The files are therefore always written as regular files. If the output path of the
    simulation includes an HDF5 file, the probe also stores each completed file as a text column
    file bundle at the end of the run, leaving the regular file in place. Only the root process
    writes output; the series hold identical values on all processes. */
class HistoryProbe : public Probe
{
    ITEM_CONCRETE(HistoryProbe, Probe, "convergence: iteration history")
        ATTRIBUTE_TYPE_DISPLAYED_IF(HistoryProbe, "IteratePrimary|IterateSecondary")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** The destructor closes the file of the current loop, if any. */
    ~HistoryProbe();

protected:
    /** This function locates the iteration history. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

protected:
    /** This function returns \c Iterations, indicating that the probe is performed after each
        iteration of each iteration loop. */
    When when() const override;

    /** This function opens the output file for the current loop when it is first performed in that
        loop, and appends a row for the current iteration. */
    void probe() override;

public:
    /** This function closes the file of the last loop, and, if the output path of the simulation
        includes an HDF5 file, stores each file written by the probe as a text column file bundle.
        It then calls the implementation in the base class. */
    void probeRun() override;

    //======================== Private helpers =======================

private:
    /** This function closes the file of the current loop, if any. */
    void closeFile();

    //======================== Data Members ========================

private:
    IterationHistory* _history{nullptr};
    std::unique_ptr<TextOutFile> _file;                              // the file for the current loop, or null
    IterationHistory::Loop _fileLoop{IterationHistory::Loop::None};  // the loop for the current file
    vector<const ScalarSeries*> _columns;  // the series written to the current file, in column order
    vector<string> _fileNames;             // the names of the files written so far, without prefix and extension
};

////////////////////////////////////////////////////////////////////

#endif
