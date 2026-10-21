/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "HistoryProbe.hpp"
#include "ColumnOutFile.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "TextOutFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the name of the specified iteration loop, for use in file names and descriptions
    string loopName(IterationHistory::Loop loop)
    {
        switch (loop)
        {
            case IterationHistory::Loop::Primary: return "primary";
            case IterationHistory::Loop::Secondary: return "secondary";
            case IterationHistory::Loop::Merged: return "merged";
            case IterationHistory::Loop::None: break;
        }
        return "none";
    }
}

////////////////////////////////////////////////////////////////////

HistoryProbe::~HistoryProbe()
{
    closeFile();
}

////////////////////////////////////////////////////////////////////

void HistoryProbe::setupSelfBefore()
{
    Probe::setupSelfBefore();
    _history = find<IterationHistory>();
}

////////////////////////////////////////////////////////////////////

Probe::When HistoryProbe::when() const
{
    return When::Iterations;
}

////////////////////////////////////////////////////////////////////

void HistoryProbe::probe()
{
    auto units = find<Units>();
    auto loop = _history->loop();

    // open the file for this loop when the probe is first performed in the loop
    if (!_file || _fileLoop != loop)
    {
        closeFile();
        _fileLoop = loop;
        _columns = _history->allScalarSeries();

        // locate the media, if any, to describe the series kept on behalf of a medium component by its index
        vector<Medium*> media;
        if (auto ms = find<MediumSystem>(false)) media = ms->media();

        string fileName = probeName() + "_" + loopName(loop);
        _fileNames.push_back(fileName);
        _file = std::make_unique<TextOutFile>(this, fileName, loopName(loop) + " iteration history");
        _file->writeLine("# Convergence history of the " + loopName(loop) + " iteration loop");
        _file->addColumn("iteration index", "", 'd');
        for (auto series : _columns)
        {
            auto item = series->key().item;
            auto position = std::find(media.begin(), media.end(), item);
            string itemDescription = position != media.end() ? "medium " + std::to_string(position - media.begin())
                                                             : (item ? item->type() : string("unknown item"));
            string unit = series->quantity().empty() ? "" : units->unit(series->quantity());
            _file->addColumn(series->description() + " -- " + itemDescription, unit, 'g');
        }
    }

    // warn about series that were declared after the file has been opened
    auto all = _history->allScalarSeries();
    if (all.size() > _columns.size())
        find<Log>()->warning(typeAndName() + " leaves out " + std::to_string(all.size() - _columns.size())
                             + " series declared after the first iteration of the " + loopName(loop) + " loop");

    // append a row for the current iteration, and flush the file
    vector<double> row;
    row.push_back(_history->iteration());
    for (auto series : _columns)
    {
        double value = series->value();
        row.push_back(series->quantity().empty() ? value : units->out(series->quantity(), value));
    }
    _file->writeRow(row);
    _file->flush();
}

////////////////////////////////////////////////////////////////////

void HistoryProbe::probeRun()
{
    // complete the file of the last loop, and store all files as bundles if the output path includes an HDF5 file
    closeFile();
    for (const auto& fileName : _fileNames) ColumnOutFile::copyToBundle(this, fileName, "iteration history");

    Probe::probeRun();
}

////////////////////////////////////////////////////////////////////

void HistoryProbe::closeFile()
{
    if (_file)
    {
        _file->close();
        _file.reset();
    }
}

////////////////////////////////////////////////////////////////////
