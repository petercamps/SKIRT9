/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EscapeInstrument.hpp"
#include "ColumnOutFile.hpp"
#include "FluxRecorder.hpp"
#include "PhotonPacket.hpp"

////////////////////////////////////////////////////////////////////

void EscapeInstrument::setupSelfBefore()
{
    Instrument::setupSelfBefore();

    // configure the flux recorder with the appropriate frame and distances
    configureRecorderDistance(distance());
}

////////////////////////////////////////////////////////////////////

void EscapeInstrument::setupSelfAfter()
{
    // configure the flux recorder for the escape bins; the binning has been set up at this point
    int numBins = _binning->numBins();
    vector<double> skyFractions(numBins);
    vector<string> descriptions(numBins);
    for (int b = 0; b != numBins; ++b)
    {
        skyFractions[b] = _binning->skyFraction(b);
        descriptions[b] = _binning->describeBin(b);
    }
    instrumentFluxRecorder()->includeFluxDensityForEscapeBins(skyFractions, descriptions);

    // for the data cube style, arrange the bins in a frame as provided by the binning
    if (_outputStyle == OutputStyle::DataCubes)
    {
        auto layout = _binning->frameLayout();
        vector<int> pixelIndices(numBins);
        for (int b = 0; b != numBins; ++b) pixelIndices[b] = _binning->framePixel(b);
        instrumentFluxRecorder()->includeDataCubesForEscapeBins(layout.numPixelsX, layout.numPixelsY, pixelIndices,
                                                                layout.incrementX, layout.incrementY, layout.centerX,
                                                                layout.centerY, layout.quantityXY);
    }

    // finalize the configuration of the flux recorder
    Instrument::setupSelfAfter();
}

////////////////////////////////////////////////////////////////////

void EscapeInstrument::detectEscape(const PhotonPacket* pp, double escapeFraction, bool includeTransparent)
{
    if (!recordsWavelength(pp->wavelength())) return;
    int bin = _binning->bin(pp->direction());
    if (bin < 0) return;
    instrumentFluxRecorder()->detect(pp, bin, std::numeric_limits<double>::infinity(), escapeFraction,
                                     includeTransparent);
}

////////////////////////////////////////////////////////////////////

void EscapeInstrument::write()
{
    Instrument::write();

    // write the bin table, with the columns provided by the binning
    ColumnOutFile file(this, instrumentName() + "_bins", "escape bins");
    file.addColumn("bin index", "", 'd');
    file.addColumn("sky fraction");
    _binning->addBinTableColumns(file);
    int numBins = _binning->numBins();
    for (int b = 0; b != numBins; ++b)
    {
        vector<double> values({static_cast<double>(b), _binning->skyFraction(b)});
        for (double value : _binning->binTableValues(b)) values.push_back(value);
        file.writeRow(values);
    }
    file.close();
}

////////////////////////////////////////////////////////////////////
