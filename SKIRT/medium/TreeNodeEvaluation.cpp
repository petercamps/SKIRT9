/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TreeNodeEvaluation.hpp"
#include "MassInBoxInterface.hpp"
#include "Medium.hpp"
#include "MediumSystem.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

TreeNodeEvaluation::TreeNodeEvaluation(const MediumSystem* ms, Random* random, int numSamples)
    : _ms(ms), _random(random), _numSamples(numSamples)
{
    if (_ms)
    {
        _media = _ms->media();
        for (auto medium : _media)
        {
            _mibv.push_back(medium->interface<MassInBoxInterface>(0, 0, false));
            _dustv.push_back(medium->mix()->materialType() == MaterialMix::MaterialType::Dust);
        }
    }
    int numMedia = _media.size();
    _hasSamplesv.resize(numMedia);
    _samplesv.resize(numMedia);
    _hasAmountv.resize(numMedia);
    _amountv.resize(numMedia);
}

////////////////////////////////////////////////////////////////////

void TreeNodeEvaluation::reset(const Box& box, int level)
{
    _box = box;
    _level = level;
    _hasPositions = false;
    std::fill(_hasSamplesv.begin(), _hasSamplesv.end(), 0);
    std::fill(_hasAmountv.begin(), _hasAmountv.end(), 0);
    _hasTypeSamplesv.fill(false);
}

////////////////////////////////////////////////////////////////////

double TreeNodeEvaluation::density(int h)
{
    if (_mibv[h]) return componentAmount(h) / _box.volume();

    const auto& samples = componentDensitySamples(h);
    return std::accumulate(samples.cbegin(), samples.cend(), 0.) / _numSamples;
}

////////////////////////////////////////////////////////////////////

double TreeNodeEvaluation::amount(MaterialMix::MaterialType type)
{
    double amount = 0.;
    for (int h : indices(type)) amount += componentAmount(h);
    return amount;
}

////////////////////////////////////////////////////////////////////

const vector<double>& TreeNodeEvaluation::densitySamples(MaterialMix::MaterialType type)
{
    // for a single medium component, simply return its samples
    const auto& hv = indices(type);
    if (hv.size() == 1) return componentDensitySamples(hv[0]);

    // otherwise, sum the samples of the medium components
    int t = static_cast<int>(type);
    auto& samples = _typeSamplesv[t];
    if (!_hasTypeSamplesv[t])
    {
        samples.assign(_numSamples, 0.);
        for (int h : hv)
        {
            const auto& componentSamplesv = componentDensitySamples(h);
            for (int i = 0; i != _numSamples; ++i) samples[i] += componentSamplesv[i];
        }
        _hasTypeSamplesv[t] = true;
    }
    return samples;
}

////////////////////////////////////////////////////////////////////

string TreeNodeEvaluation::materialTypeName(MaterialMix::MaterialType type)
{
    switch (type)
    {
        case MaterialMix::MaterialType::Dust: return "dust";
        case MaterialMix::MaterialType::Electrons: return "electrons";
        case MaterialMix::MaterialType::Gas: return "gas";
    }
    return string();
}

////////////////////////////////////////////////////////////////////

const vector<int>& TreeNodeEvaluation::indices(MaterialMix::MaterialType type) const
{
    return _ms ? _ms->mediumIndices(type) : _noIndices;
}

////////////////////////////////////////////////////////////////////

const vector<Position>& TreeNodeEvaluation::samplePositions()
{
    if (!_hasPositions)
    {
        _positionv.resize(_numSamples);
        for (auto& position : _positionv) position = _random->position(_box);
        _hasPositions = true;
    }
    return _positionv;
}

////////////////////////////////////////////////////////////////////

const vector<double>& TreeNodeEvaluation::componentDensitySamples(int h)
{
    auto& samples = _samplesv[h];
    if (!_hasSamplesv[h])
    {
        const auto& positionv = samplePositions();
        samples.resize(_numSamples);
        for (int i = 0; i != _numSamples; ++i)
            samples[i] = _dustv[h] ? _media[h]->massDensity(positionv[i]) : _media[h]->numberDensity(positionv[i]);
        _hasSamplesv[h] = true;
    }
    return samples;
}

////////////////////////////////////////////////////////////////////

double TreeNodeEvaluation::componentAmount(int h)
{
    if (!_hasAmountv[h])
    {
        if (_mibv[h])
        {
            _amountv[h] = _dustv[h] ? _mibv[h]->massInBox(_box) : _mibv[h]->numberInBox(_box);
        }
        else
        {
            const auto& samples = componentDensitySamples(h);
            _amountv[h] = std::accumulate(samples.cbegin(), samples.cend(), 0.) / _numSamples * _box.volume();
        }
        _hasAmountv[h] = true;
    }
    return _amountv[h];
}

////////////////////////////////////////////////////////////////////
