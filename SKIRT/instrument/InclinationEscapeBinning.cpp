/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "InclinationEscapeBinning.hpp"
#include "ColumnOutFile.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

int InclinationEscapeBinning::numBins() const
{
    return _numInclinationBins * _numAzimuthBins;
}

////////////////////////////////////////////////////////////////////

int InclinationEscapeBinning::bin(Direction bfk) const
{
    // the inclination bins are equally spaced in cos(theta), or in |cos(theta)| if the hemispheres are folded
    double fraction = _foldHemispheres ? 1. - abs(bfk.kz()) : 0.5 * (1. - bfk.kz());
    int i = std::clamp(static_cast<int>(fraction * _numInclinationBins), 0, _numInclinationBins - 1);

    // the azimuth bins have equal width in phi, starting at phi = 0
    int k = 0;
    if (_numAzimuthBins > 1)
    {
        double phi = atan2(bfk.ky(), bfk.kx());
        if (phi < 0.) phi += 2. * M_PI;
        k = std::clamp(static_cast<int>(phi / (2. * M_PI) * _numAzimuthBins), 0, _numAzimuthBins - 1);
    }
    return i * _numAzimuthBins + k;
}

////////////////////////////////////////////////////////////////////

double InclinationEscapeBinning::skyFraction(int /*bin*/) const
{
    return 1. / numBins();
}

////////////////////////////////////////////////////////////////////

std::tuple<double, double, double, double> InclinationEscapeBinning::binRanges(int bin) const
{
    int i = bin / _numAzimuthBins;
    int k = bin % _numAzimuthBins;
    double width = _foldHemispheres ? 1. / _numInclinationBins : 2. / _numInclinationBins;
    return {acos(1. - i * width), acos(max(-1., 1. - (i + 1) * width)), 2. * M_PI * k / _numAzimuthBins,
            2. * M_PI * (k + 1) / _numAzimuthBins};
}

////////////////////////////////////////////////////////////////////

string InclinationEscapeBinning::describeBin(int bin) const
{
    auto units = find<Units>();
    auto angle = [units](double value) {
        return StringUtils::toString(units->oposangle(value), 'g', 4) + " " + units->uposangle();
    };

    auto [thetaMin, thetaMax, phiMin, phiMax] = binRanges(bin);
    string description = "inclination " + angle(thetaMin) + " to " + angle(thetaMax);
    if (_foldHemispheres) description += " in both hemispheres";
    if (_numAzimuthBins > 1) description += ", azimuth " + angle(phiMin) + " to " + angle(phiMax);
    return description;
}

////////////////////////////////////////////////////////////////////

void InclinationEscapeBinning::addBinTableColumns(ColumnOutFile& file) const
{
    auto units = find<Units>();
    file.addColumn("minimum inclination", units->uposangle());
    file.addColumn("maximum inclination", units->uposangle());
    if (_numAzimuthBins > 1)
    {
        file.addColumn("minimum azimuth", units->uposangle());
        file.addColumn("maximum azimuth", units->uposangle());
    }
}

////////////////////////////////////////////////////////////////////

vector<double> InclinationEscapeBinning::binTableValues(int bin) const
{
    auto units = find<Units>();
    auto [thetaMin, thetaMax, phiMin, phiMax] = binRanges(bin);
    vector<double> values({units->oposangle(thetaMin), units->oposangle(thetaMax)});
    if (_numAzimuthBins > 1)
    {
        values.push_back(units->oposangle(phiMin));
        values.push_back(units->oposangle(phiMax));
    }
    return values;
}

////////////////////////////////////////////////////////////////////
