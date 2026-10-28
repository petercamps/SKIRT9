/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumSystem.hpp"
#include "Configuration.hpp"
#include "DensityInCellInterface.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "HistoryScope.hpp"
#include "IterationHistory.hpp"
#include "LockFree.hpp"
#include "Log.hpp"
#include "MaterialMix.hpp"
#include "MaterialState.hpp"
#include "NR.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PathSegmentGenerator.hpp"
#include "PhotonPacket.hpp"
#include "ProcessManager.hpp"
#include "Random.hpp"
#include "ShortArray.hpp"
#include "SpecialFunctions.hpp"
#include "StateVariableList.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "Table.hpp"
#include "TreeSpatialGrid.hpp"
#include <unordered_set>

////////////////////////////////////////////////////////////////////

namespace
{
    // maximum number of cell densities calculated between two invocations of infoIfElapsed()
    const size_t logProgressChunkSize = 10000;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This helper class performs spatial sampling of the medium properties within a given cell and
    // for a given medium component. The number of samples to be taken can be specified separately
    // for the density and for all other properties (velocity, magnetic field, imported parameters).
    // The main reason for bundling this functionality in its own class/object is to allow caching
    // the random sampling positions with the corresponding sampled density values in a cell.
    // The constructor takes arguments that depend only on the simulation configuration and thus
    // do not vary between cells. This allows consecutively reusing the same object for multiple
    // cells, avoiding memory allocation/deallocation for each cell.
    class PropertySampler
    {
    private:
        vector<Medium*>& _media;
        SpatialGrid* _grid;
        int _numDensitySamples;
        int _numPropertySamples;
        int _numMedia;
        int _cellIndex;
        Position _center;
        vector<Position> _positions;              // indexed on n
        Table<2> _densities;                      // indexed on h and n
        vector<vector<int>> _categoricalIndices;  // indices of the categorical parameters, indexed on h

    public:
        // The constructor allocates the required memory depending on the sampling configuration.
        // The number of samples for each type (density or other properties) can be:
        //   0: the sampling function(s) will never be called
        //   1: "sample" at the cell center only
        //  >1: sample at the specified nr of random points in the cell
        PropertySampler(vector<Medium*>& media, SpatialGrid* grid, int numDensitySamples, int numPropertySamples)
            : _media(media), _grid(grid), _numDensitySamples(numDensitySamples),
              _numPropertySamples(numPropertySamples), _numMedia(_media.size()), _cellIndex(0)
        {
            if (_numDensitySamples > 1 || _numPropertySamples > 1)
            {
                int numSamples = max(_numDensitySamples, _numPropertySamples);
                _positions.resize(numSamples);
                _densities.resize(_numMedia, numSamples);
            }

            // if we need random sampling of the custom parameters, remember which ones are categorical
            if (_numPropertySamples > 1)
            {
                _categoricalIndices.resize(_numMedia);
                for (int h = 0; h != _numMedia; ++h)
                {
                    auto info = _media[h]->parameterInfo();
                    int numParams = info.size();
                    for (int i = 0; i != numParams; ++i)
                        if (info[i].kind() == QuantityKind::Categorical) _categoricalIndices[h].push_back(i);
                }
            }
        }

        // This function should be called with a spatial cell index before calling any of the other
        // functions for that cell. It initializes the required sampling positions and obtains the
        // corresponding sampled density values.
        void prepareForCell(int cellIndex)
        {
            _cellIndex = cellIndex;

            // if we need center "sampling", get the cell center
            if (_numDensitySamples == 1 || _numPropertySamples == 1)
            {
                _center = _grid->centralPositionInCell(_cellIndex);
            }

            // if we need random sampling, draw the positions and get the corresponding densities
            int numSamples = _positions.size();
            for (int n = 0; n != numSamples; ++n)
            {
                _positions[n] = _grid->randomPositionInCell(_cellIndex);
                for (int h = 0; h != _numMedia; ++h) _densities(h, n) = _media[h]->numberDensity(_positions[n]);
            }
        }

        // Returns the central or average density in the cell for the given medium component.
        double density(int h)
        {
            if (_numDensitySamples == 1) return _media[h]->numberDensity(_center);
            double sum = 0.;
            for (int n = 0; n != _numDensitySamples; ++n) sum += _densities(h, n);
            return sum / _numDensitySamples;
        }

        // Returns the central or density-averaged velocity in the cell for the given medium component.
        Vec bulkVelocity(int h)
        {
            if (_numPropertySamples == 1) return _media[h]->bulkVelocity(_center);
            double nsum = 0;
            Vec vsum;
            for (int n = 0; n != _numPropertySamples; ++n)
            {
                nsum += _densities(h, n);
                vsum += _densities(h, n) * _media[h]->bulkVelocity(_positions[n]);
            }
            return nsum > 0. ? vsum / nsum : Vec();
        }

        // Returns the central or average magnetic field in the cell for the given medium component.
        // Note the magnetic field average is NOT weighed with density.
        Vec magneticField(int h)
        {
            if (_numPropertySamples == 1) return _media[h]->magneticField(_center);
            Vec Bsum;
            for (int n = 0; n != _numPropertySamples; ++n) Bsum += _media[h]->magneticField(_positions[n]);
            return Bsum / _numPropertySamples;
        }

        // Returns the central or sampled custom parameter values in the cell for the given medium component:
        // the density-averaged value for an intensive parameter, and the value of the sample with the highest
        // density (the first such sample if several samples have this density) for a categorical parameter.
        void parameters(int h, Array& params)
        {
            // always sample at the center to make sure that the output array has the proper size
            _media[h]->parameters(_center, params);
            // if center-sampling is requested, we're done
            if (_numPropertySamples == 1) return;
            // otherwise, clear the output values and accumulate the samples, remembering the densest sample
            params = 0.;
            double nsum = 0;
            double nmax = -1.;
            Array sample, densest;
            for (int n = 0; n != _numPropertySamples; ++n)
            {
                nsum += _densities(h, n);
                _media[h]->parameters(_positions[n], sample);
                params += _densities(h, n) * sample;
                if (!_categoricalIndices[h].empty() && _densities(h, n) > nmax)
                {
                    nmax = _densities(h, n);
                    densest = sample;
                }
            }
            if (nsum > 0.)
                params /= nsum;
            else
                params = 0.;
            for (int i : _categoricalIndices[h]) params[i] = densest[i];
            return;
        }
    };
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the physical quantity of the volume integral of a state variable with the specified quantity, as a
    // quantity name known to the Units class, or the empty string for a dimensionless or unknown quantity
    string integratedQuantity(const string& quantity)
    {
        if (quantity == "massvolumedensity") return "mass";
        return string();  // the volume integral of a number density is a pure number
    }
}

////////////////////////////////////////////////////////////////////

void MediumSystem::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    // set up the children that precede the spatial grid, in the usual order (some medium components need the
    // options to be set up first); the remaining children (i.e. the spatial grid) are set up after this function
    for (Item* child : children())
    {
        if (child == _grid) break;
        if (auto item = dynamic_cast<SimulationItem*>(child)) item->setup();
    }

    // cache a list of medium component indices for each material type, so that these are available to the grid
    int numMedia = _media.size();
    for (int h = 0; h != numMedia; ++h)
    {
        switch (_media[h]->capabilities().materialType)
        {
            case MaterialMix::MaterialType::Dust: _dust_hv.push_back(h); break;
            case MaterialMix::MaterialType::Gas: _gas_hv.push_back(h); break;
            case MaterialMix::MaterialType::Electrons: _elec_hv.push_back(h); break;
        }
    }
}

////////////////////////////////////////////////////////////////////

void MediumSystem::setupSelfAfter()
{
    SimulationItem::setupSelfAfter();
    auto log = find<Log>();
    auto parfac = find<ParallelFactory>();
    _config = find<Configuration>();
    _history = find<IterationHistory>();

    _numCells = _grid->numCells();
    if (_numCells < 1) throw FATALERROR("The spatial grid must have at least one cell");
    _numMedia = _media.size();
    size_t allocatedBytes = 0;

    // ----- allocate memory for the medium state -----

    // basic configuration
    _state.initConfiguration(_numCells, _numMedia);

    // common state variables
    vector<StateVariable> variables;
    variables.emplace_back(StateVariable::volume());
    if (_config->hasMovingMedia()) variables.emplace_back(StateVariable::bulkVelocity());
    if (_config->hasMagneticField()) variables.emplace_back(StateVariable::magneticField());
    _state.initCommonStateVariables(variables);

    // specific state variables; for each variable that the material mix marked for aggregation, declare an aggregate
    // series in the iteration history on behalf of the mix, keyed on the medium component and the id chosen by the mix
    // (for an array variable, a series for each element, with consecutive ids)
    for (auto medium : _media)
    {
        const auto& variables = medium->mix()->specificStateVariableInfo();
        _state.initSpecificStateVariables(variables);
        for (const StateVariable& variable : variables)
        {
            if (variable.isAggregated())
            {
                for (int k = 0; k != variable.count(); ++k)
                    _history->aggregateSeries(medium, variable.aggregateSeriesId() + k,
                                              AggregationRule{medium, variable.customIndex() + k}, 2,
                                              variable.elementDescription(k), integratedQuantity(variable.quantity()));
            }
        }
    }

    // finalize
    allocatedBytes += _state.initAllocate() * sizeof(double);

    // ----- allocate memory for the radiation field -----

    if (_config->hasRadiationField())
    {
        _wavelengthGrid = _config->radiationFieldWLG();
        _rf1.resize(_numCells, _wavelengthGrid->numBins());
        allocatedBytes += _rf1.size() * sizeof(double);

        if (_config->hasSecondaryRadiationField())
        {
            _rf2.resize(_numCells, _wavelengthGrid->numBins());
            _rf2c.resize(_numCells, _wavelengthGrid->numBins());
            allocatedBytes += 2 * _rf2.size() * sizeof(double);
        }

        // allocate memory for the dust-absorbed luminosity per cell, following the corresponding tables;
        // use the list of dust components rather than hasDust(), because the material mixes are obtained below
        if (_config->hasDustEmission() && !_dust_hv.empty())
        {
            _labs1.resize(_numCells);
            allocatedBytes += _labs1.size() * sizeof(double);

            if (_config->hasSecondaryRadiationField())
            {
                _labs2.resize(_numCells);
                _labs2c.resize(_numCells);
                allocatedBytes += 2 * _labs2.size() * sizeof(double);
            }
        }
    }

    // ----- create the self-heating rebalancer -----

    if (_config->rebalanceSecondaryEmission() && _labs2.size())
    {
        _rebalancer = std::make_unique<SelfHeatingRebalancer>(
            _config->numAbsorptionRatioBins(), _config->numHeatingRateBins(), _history, dustEmissionOptions(), log);
    }

    // ----- cache info on the dust emission wavelength grid -----

    if (_config->hasDustEmission())
    {
        _numDustEmissionWavelengths = _config->dustEmissionWLG()->extlambdav().size();
    }

    // ----- obtain the material mix pointers -----

    if (_config->hasVariableMedia())
    {
        // we need a separate material mix pointer for each spatial cell (per medium component)
        _mixv.resize(_numCells * _numMedia);
        _mixPerCell = true;
        for (int m = 0; m != _numCells; ++m)
        {
            Position bfr = _grid->centralPositionInCell(m);
            for (int h = 0; h != _numMedia; ++h) _mixv[m * _numMedia + h] = _media[h]->mix(bfr);
        }

        // verify that each distinct material mix of a medium component with a variable mix has the capabilities
        // of the medium component and the specific state variables of its representative mix, which determine
        // the configuration and the layout of the medium state, and, if the medium component declares a uniform
        // mass per entity, the mass of its representative mix, which converts between masses and numbers
        for (int h = 0; h != _numMedia; ++h)
        {
            if (!_media[h]->hasVariableMix()) continue;
            auto capabilities = _media[h]->capabilities();
            const auto& variables = _media[h]->mix()->specificStateVariableInfo();
            bool uniformMass = _media[h]->hasUniformMass();
            double mass = _media[h]->mix()->mass();
            std::unordered_set<const MaterialMix*> verified;
            for (int m = 0; m != _numCells; ++m)
            {
                const MaterialMix* mix = _mixv[m * _numMedia + h];
                if (!verified.insert(mix).second) continue;
                if (mix->capabilities() != capabilities)
                    throw FATALERROR("A material mix of medium component " + std::to_string(h)
                                     + " has other capabilities than the medium component");
                if (mix->specificStateVariableInfo() != variables)
                    throw FATALERROR("A material mix of medium component " + std::to_string(h)
                                     + " has other specific state variables than the representative mix");
                if (uniformMass && mix->mass() != mass)
                    throw FATALERROR("A material mix of medium component " + std::to_string(h)
                                     + " has another mass per entity than the representative mix, although its"
                                     + " material mix family declares a uniform mass");
            }
        }
    }
    else
    {
        // the material mix pointer is identical for all spatial cells (per medium component)
        _mixv.resize(_numMedia);
        for (int h = 0; h != _numMedia; ++h) _mixv[h] = _media[h]->mix();
    }
    allocatedBytes += _mixv.size() * sizeof(MaterialMix*);

    // cache a list of medium component indices for each type of dynamic medium state
    for (int h = 0; h != _numMedia; ++h)
    {
        switch (mix(0, h)->hasDynamicMediumState())
        {
            case MaterialMix::DynamicStateType::None: break;
            case MaterialMix::DynamicStateType::Primary: _pdms_hv.push_back(h); break;
            case MaterialMix::DynamicStateType::Secondary: _sdms_hv.push_back(h); break;
            case MaterialMix::DynamicStateType::PrimaryIfMergedIterations:
                if (_config->hasMergedIterations())
                    _pdms_hv.push_back(h);
                else
                    _sdms_hv.push_back(h);
                break;
        }
    }

    // ----- inform user about allocated memory -----

    log->info(typeAndName() + " allocated " + StringUtils::toMemSizeString(allocatedBytes) + " of memory");

    // ----- calculate medium properties parallelized on spatial cells -----

    log->info("Determining medium properties for " + std::to_string(_numCells) + " cells...");
    auto dic = _grid->interface<DensityInCellInterface>(0, 0, false);  // optional fast-track interface for densities
    if (dic) log->info("  (obtaining densities through calculation rather than sampling)");
    log->infoSetElapsed(_numCells);
    parfac->parallelDistributed()->call(_numCells, [this, log, dic](size_t firstIndex, size_t numIndices) {
        // construct a property sampler to be shared by all cells handled in the loop below
        bool hasProperty = _config->hasMovingMedia() || _config->hasMagneticField();
        for (int h = 0; h != _numMedia; ++h)
            if (_media[h]->hasParameters()) hasProperty = true;
        PropertySampler sampler(_media, _grid, dic ? 0 : _config->numDensitySamples(),
                                hasProperty ? _config->numPropertySamples() : 0);

        // loop over cells
        while (numIndices)
        {
            size_t currentChunkSize = min(logProgressChunkSize, numIndices);
            for (size_t m = firstIndex; m != firstIndex + currentChunkSize; ++m)
            {
                // prepare the sampler for this cell (i.e. store the relevant positions and density samples)
                sampler.prepareForCell(m);

                // volume
                _state.setVolume(m, _grid->volume(m));

                // density: use optional fast-track interface or sample within the cell
                if (dic)
                {
                    for (int h = 0; h != _numMedia; ++h) _state.setNumberDensity(m, h, dic->numberDensity(h, m));
                }
                else
                {
                    for (int h = 0; h != _numMedia; ++h) _state.setNumberDensity(m, h, sampler.density(h));
                }

                // bulk velocity: aggregate values sampled for each component according to specfied policy
                if (_config->hasMovingMedia())
                {
                    // accumulate total density in the cell
                    double n = 0.;
                    for (int h = 0; h != _numMedia; ++h) n += _state.numberDensity(m, h);

                    // leave bulk velocity at zero if cell has no material
                    if (n > 0.)
                    {
                        Vec v;
                        switch (_samplingOptions->aggregateVelocity())
                        {
                            case SamplingOptions::AggregatePolicy::Average:
                                for (int h = 0; h != _numMedia; ++h)
                                {
                                    v += _state.numberDensity(m, h) * sampler.bulkVelocity(h);
                                }
                                v /= n;
                                break;
                            case SamplingOptions::AggregatePolicy::Maximum:
                                for (int h = 0; h != _numMedia; ++h)
                                {
                                    if (_state.numberDensity(m, h) > 0)
                                    {
                                        Vec w = sampler.bulkVelocity(h);
                                        if (w.norm2() > v.norm2()) v = w;
                                    }
                                }
                                break;
                            case SamplingOptions::AggregatePolicy::First:
                                for (int h = 0; h != _numMedia; ++h)
                                {
                                    if (_media[h]->hasVelocity())
                                    {
                                        if (_state.numberDensity(m, h) > 0) v = sampler.bulkVelocity(h);
                                        break;
                                    }
                                    break;
                                }
                        }
                        _state.setBulkVelocity(m, v);
                    }
                }

                // magnetic field: retrieve value sampled from medium component that specifies it
                if (_config->hasMagneticField())
                {
                    _state.setMagneticField(m, sampler.magneticField(_config->magneticFieldMediumIndex()));
                }

                // specific state variables other than density:
                // retrieve the imported parameters sampled from the corresponding medium component, except for a
                // medium with a variable mix, for which the parameters have already selected the material mix
                for (int h = 0; h != _numMedia; ++h)
                {
                    Array params;
                    if (_media[h]->hasParameters() && !_media[h]->hasVariableMix()) sampler.parameters(h, params);
                    MaterialState mst(_state, m, h);
                    mix(m, h)->initializeSpecificState(&mst, params);
                }
            }
            log->infoIfElapsed("Determined medium properties: ", currentChunkSize);
            firstIndex += currentChunkSize;
            numIndices -= currentChunkSize;
        }
    });

    // communicate the calculated states across multiple processes, if needed
    _state.initCommunicate();

    // record the aggregates of the initial state
    recordAggregates();

    log->info("Done determining medium properties");

    // verify the dynamic refinement configuration and reserve memory, if applicable
    setupDynamicRefinement();
}

////////////////////////////////////////////////////////////////////

int MediumSystem::dimension() const
{
    int result = 1;
    for (auto medium : _media) result = max(result, medium->dimension());
    return result;
}

////////////////////////////////////////////////////////////////////

int MediumSystem::gridDimension() const
{
    return _grid->dimension();
}

////////////////////////////////////////////////////////////////////

int MediumSystem::numMedia() const
{
    return _numMedia;
}

////////////////////////////////////////////////////////////////////

int MediumSystem::numCells() const
{
    return _numCells;
}

////////////////////////////////////////////////////////////////////

const MaterialMix* MediumSystem::mix(int m, int h) const
{
    return _mixPerCell ? _mixv[m * _numMedia + h] : _mixv[h];
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::hasMaterialType(MaterialMix::MaterialType type) const
{
    for (int h = 0; h != _numMedia; ++h)
        if (mix(0, h)->materialType() == type) return true;
    return false;
}

////////////////////////////////////////////////////////////////////

const vector<int>& MediumSystem::mediumIndices(MaterialMix::MaterialType type) const
{
    switch (type)
    {
        case MaterialMix::MaterialType::Dust: return _dust_hv;
        case MaterialMix::MaterialType::Electrons: return _elec_hv;
        case MaterialMix::MaterialType::Gas: return _gas_hv;
    }
    return _dust_hv;  // unreachable; avoids a compiler warning
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::isMaterialType(MaterialMix::MaterialType type, int h) const
{
    return mix(0, h)->materialType() == type;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::dustMassDensity(Position bfr) const
{
    double result = 0.;
    for (int h : _dust_hv) result += _media[h]->massDensity(bfr);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::electronNumberDensity(Position bfr) const
{
    double result = 0.;
    for (int h : _elec_hv) result += _media[h]->numberDensity(bfr);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::gasNumberDensity(Position bfr) const
{
    double result = 0.;
    for (int h : _gas_hv) result += _media[h]->numberDensity(bfr);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::volume(int m) const
{
    return _state.volume(m);
}

////////////////////////////////////////////////////////////////////

Vec MediumSystem::bulkVelocity(int m) const
{
    return _state.bulkVelocity(m);
}

////////////////////////////////////////////////////////////////////

Vec MediumSystem::magneticField(int m) const
{
    return _state.magneticField(m);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::massDensity(int m) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h) result += massDensity(m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::dustMassDensity(int m) const
{
    double result = 0.;
    for (int h : _dust_hv) result += massDensity(m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::electronNumberDensity(int m) const
{
    double result = 0.;
    for (int h : _elec_hv) result += numberDensity(m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::gasNumberDensity(int m) const
{
    double result = 0.;
    for (int h : _gas_hv) result += numberDensity(m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::numberDensity(int m, int h) const
{
    return _state.numberDensity(m, h);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::massDensity(int m, int h) const
{
    return _state.numberDensity(m, h) * mix(m, h)->mass();
}

////////////////////////////////////////////////////////////////////

double MediumSystem::custom(int m, int h, int i) const
{
    return _state.custom(m, h, i);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, int h) const
{
    MaterialState mst(_state, m, h);
    return mix(m, h)->opacityAbs(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacitySca(double lambda, int m, int h) const
{
    MaterialState mst(_state, m, h);
    return mix(m, h)->opacitySca(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, int h) const
{
    MaterialState mst(_state, m, h);
    return mix(m, h)->opacityExt(lambda, &mst, nullptr);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, MaterialMix::MaterialType type) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
        if (mix(0, h)->materialType() == type) result += opacityAbs(lambda, m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, MaterialMix::MaterialType type) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
        if (mix(0, h)->materialType() == type) result += opacityExt(lambda, m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h) result += opacityExt(lambda, m, h);
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityAbs(double lambda, int m, const PhotonPacket* pp) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        result += mix(m, h)->opacityAbs(lambda, &mst, pp);
    }
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacitySca(double lambda, int m, const PhotonPacket* pp) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        result += mix(m, h)->opacitySca(lambda, &mst, pp);
    }
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::opacityExt(double lambda, int m, const PhotonPacket* pp) const
{
    double result = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        result += mix(m, h)->opacityExt(lambda, &mst, pp);
    }
    return result;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::perceivedWavelengthForScattering(const PhotonPacket* pp) const
{
    if (_config->hasMovingMedia())
        return pp->perceivedWavelength(_state.bulkVelocity(pp->interactionCellIndex()),
                                       _config->hubbleExpansionRate() * pp->interactionDistance());
    else
        return pp->wavelength();
}

////////////////////////////////////////////////////////////////////

double MediumSystem::albedoForScattering(const PhotonPacket* pp) const
{
    double lambda = perceivedWavelengthForScattering(pp);
    int m = pp->interactionCellIndex();
    if (m < 0) throw FATALERROR("Cannot locate photon packet interaction point");

    double ksca = 0.;
    double kext = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        ksca += mix(m, h)->opacitySca(lambda, &mst, pp);
        kext += mix(m, h)->opacityExt(lambda, &mst, pp);
    }
    return kext > 0. ? ksca / kext : 0.;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::weightsForScattering(ShortArray& wv, double lambda, const PhotonPacket* pp) const
{
    // resize the target array
    wv.resize(_numMedia);

    // for a single component, the weight is trivial
    if (_numMedia == 1)
    {
        wv[0] = 1.;
        return true;
    }

    // locate the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // calculate the weights and their sum
    double sum = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        MaterialState mst(_state, m, h);
        wv[h] = mix(m, h)->opacitySca(lambda, &mst, pp);
        sum += wv[h];
    }

    // normalize the weights
    if (sum > 0.)
    {
        for (int h = 0; h != _numMedia; ++h) wv[h] /= sum;
        return true;
    }

    // none of the media scatter
    return false;
}

////////////////////////////////////////////////////////////////////

void MediumSystem::peelOffScattering(const ShortArray& wv, double lambda, Direction bfkobs, Direction bfky,
                                     PhotonPacket* pp, PhotonPacket* ppp) const
{
    // get the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // calculate the weighted sum of the effects on the Stokes vector for all media
    // we assume the wavelength does not change
    double I = 0., Q = 0., U = 0., V = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        // skip media that don't scatter this photon packet
        if (wv[h] > 0.)
        {
            double Ih = 0., Qh = 0., Uh = 0., Vh = 0.;
            MaterialState mst(_state, m, h);
            pp->setScatteringComponent(h);
            mix(m, h)->peeloffScattering(Ih, Qh, Uh, Vh, lambda, bfkobs, bfky, &mst, pp);
            I += Ih * wv[h];
            Q += Qh * wv[h];
            U += Uh * wv[h];
            V += Vh * wv[h];
        }
    }

    // pass the result to the peel-off photon packet
    ppp->launchScatteringPeelOff(pp, bfkobs, _state.bulkVelocity(m), lambda, I);
    if (_config->hasPolarization())
    {
        // the Stokes vector has already been rotated towards the observer frame,
        // so the reference direction is no longer used and can be left unspecified
        ppp->setPolarized(I, Q, U, V);
    }
}

////////////////////////////////////////////////////////////////////

void MediumSystem::peelOffScattering(int h, double w, double lambda, Direction bfkobs, Direction bfky, PhotonPacket* pp,
                                     PhotonPacket* ppp) const
{
    // get the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // calculate the effects on the Stokes vector and on the wavelength for this medium component
    double I = 0., Q = 0., U = 0., V = 0.;
    MaterialState mst(_state, m, h);
    pp->setScatteringComponent(h);
    bool emulatedSecondaryEmission = mix(m, h)->peeloffScattering(I, Q, U, V, lambda, bfkobs, bfky, &mst, pp);

    // pass the result to the peel-off photon packet
    ppp->launchScatteringPeelOff(pp, bfkobs, _state.bulkVelocity(m), lambda, I * w);
    if (_config->hasPolarization())
    {
        // the Stokes vector has already been rotated towards the observer frame,
        // so the reference direction is no longer used and can be left unspecified
        ppp->setPolarized(I, Q, U, V);
    }
    if (emulatedSecondaryEmission) ppp->setEmulatedSecondaryOrigin(h);
}

////////////////////////////////////////////////////////////////////

void MediumSystem::simulateScattering(Random* random, PhotonPacket* pp) const
{
    // locate the cell hosting the scattering event
    int m = pp->interactionCellIndex();

    // calculate the perceived wavelength in the cell
    double lambda = perceivedWavelengthForScattering(pp);

    // select a medium component within that cell
    int h = 0;
    if (_numMedia > 1)
    {
        // build the cumulative distribution corresponding to the scattering opacities
        Array Xv;
        NR::cdf(Xv, _numMedia, [this, lambda, pp, m](int h) {
            MaterialState mst(_state, m, h);
            return mix(m, h)->opacitySca(lambda, &mst, pp);
        });

        // randomly select an index
        h = NR::locateClip(Xv, random->uniform());
    }

    // actually perform the scattering event for this cell and medium component
    MaterialState mst(_state, m, h);
    pp->setScatteringComponent(h);
    mix(m, h)->performScattering(lambda, &mst, pp);
}

////////////////////////////////////////////////////////////////////

namespace
{
    // This function returns a thread-local instance of the path segment generator for the specified grid
    // that is initialized to the starting position and direction of the specified path.
    // Providing a thread-local instance avoids creating a new generator for each use.
    PathSegmentGenerator* getPathSegmentGenerator(SpatialGrid* grid, const SpatialGridPath* path)
    {
        thread_local SpatialGrid* t_grid{nullptr};
        thread_local std::unique_ptr<PathSegmentGenerator> t_generator;

        if (grid != t_grid)
        {
            t_grid = grid;
            t_generator = grid->createPathSegmentGenerator();
        }
        t_generator->start(path);
        return t_generator.get();
    }
}

////////////////////////////////////////////////////////////////////

void MediumSystem::setExtinctionOpticalDepths(PhotonPacket* pp) const
{
    // determine and store the path segments in the photon packet
    auto generator = getPathSegmentGenerator(_grid, pp);
    pp->clear();
    while (generator->next())
    {
        pp->addSegment(generator->m(), generator->ds());
    }

    // calculate the cumulative optical depth and store it in the photon packet for each path segment
    double tau = 0.;

    // single medium, spatially constant cross sections
    if (_config->hasSingleConstantSectionMedium())
    {
        double section = mix(0, 0)->sectionExt(pp->wavelength());
        for (auto& segment : pp->segments())
        {
            if (segment.m() >= 0) tau += section * _state.numberDensity(segment.m(), 0) * segment.ds();
            segment.setOpticalDepth(tau);
        }
    }

    // multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantSectionMedia())
    {
        ShortArray sectionv(_numMedia);
        for (int h = 0; h != _numMedia; ++h) sectionv[h] = mix(0, h)->sectionExt(pp->wavelength());
        for (auto& segment : pp->segments())
        {
            if (segment.m() >= 0)
                for (int h = 0; h != _numMedia; ++h)
                    tau += sectionv[h] * _state.numberDensity(segment.m(), h) * segment.ds();
            segment.setOpticalDepth(tau);
        }
    }

    // spatially variable cross sections
    else
    {
        for (auto& segment : pp->segments())
        {
            if (segment.m() >= 0)
            {
                double lambda = pp->perceivedWavelength(_state.bulkVelocity(segment.m()),
                                                        _config->hubbleExpansionRate() * segment.s());
                tau += opacityExt(lambda, segment.m(), pp) * segment.ds();
            }
            segment.setOpticalDepth(tau);
        }
    }
}

////////////////////////////////////////////////////////////////////

void MediumSystem::setScatteringAndAbsorptionOpticalDepths(PhotonPacket* pp) const
{
    // determine and store the path segments in the photon packet
    auto generator = getPathSegmentGenerator(_grid, pp);
    pp->clear();
    while (generator->next())
    {
        pp->addSegment(generator->m(), generator->ds());
    }

    // calculate the cumulative optical depths and store them in the photon packet for each path segment
    double tauSca = 0.;
    double tauAbs = 0.;

    // single medium, spatially constant cross sections
    if (_config->hasSingleConstantSectionMedium())
    {
        double sectionSca = mix(0, 0)->sectionSca(pp->wavelength());
        double sectionAbs = mix(0, 0)->sectionAbs(pp->wavelength());
        for (auto& segment : pp->segments())
        {
            if (segment.m() >= 0)
            {
                double ns = _state.numberDensity(segment.m(), 0) * segment.ds();
                tauSca += sectionSca * ns;
                tauAbs += sectionAbs * ns;
            }
            segment.setOpticalDepth(tauSca, tauAbs);
        }
    }

    // multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantSectionMedia())
    {
        ShortArray sectionScav(_numMedia);
        ShortArray sectionAbsv(_numMedia);
        for (int h = 0; h != _numMedia; ++h)
        {
            sectionScav[h] = mix(0, h)->sectionSca(pp->wavelength());
            sectionAbsv[h] = mix(0, h)->sectionAbs(pp->wavelength());
        }
        for (auto& segment : pp->segments())
        {
            if (segment.m() >= 0)
                for (int h = 0; h != _numMedia; ++h)
                {
                    double ns = _state.numberDensity(segment.m(), h) * segment.ds();
                    tauSca += sectionScav[h] * ns;
                    tauAbs += sectionAbsv[h] * ns;
                }
            segment.setOpticalDepth(tauSca, tauAbs);
        }
    }

    // spatially variable cross sections
    else
    {
        for (auto& segment : pp->segments())
        {
            if (segment.m() >= 0)
            {
                double lambda = pp->perceivedWavelength(_state.bulkVelocity(segment.m()),
                                                        _config->hubbleExpansionRate() * segment.s());
                tauSca += opacitySca(lambda, segment.m(), pp) * segment.ds();
                tauAbs += opacityAbs(lambda, segment.m(), pp) * segment.ds();
            }
            segment.setOpticalDepth(tauSca, tauAbs);
        }
    }
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::setInteractionPointUsingExtinction(PhotonPacket* pp, double tauinteract) const
{
    auto generator = getPathSegmentGenerator(_grid, pp);
    double tau = 0.;
    double s = 0.;

    // loop over the segments of the path until the interaction optical depth is reached or the path ends

    // --> single medium, spatially constant cross sections
    if (_config->hasSingleConstantSectionMedium())
    {
        double section = mix(0, 0)->sectionExt(pp->wavelength());
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tau0 = tau;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0) tau += section * _state.numberDensity(m, 0) * generator->ds();
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauinteract < tau)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauinteract, tau0, tau, s0, s));
                return true;
            }
        }
    }

    // --> multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantSectionMedia())
    {
        ShortArray sectionv(_numMedia);
        for (int h = 0; h != _numMedia; ++h) sectionv[h] = mix(0, h)->sectionExt(pp->wavelength());
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tau0 = tau;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
                for (int h = 0; h != _numMedia; ++h) tau += sectionv[h] * _state.numberDensity(m, h) * ds;
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauinteract < tau)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauinteract, tau0, tau, s0, s));
                return true;
            }
        }
    }

    // --> spatially variable cross sections
    else
    {
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tau0 = tau;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                double lambda = pp->perceivedWavelength(_state.bulkVelocity(m), _config->hubbleExpansionRate() * s);
                tau += opacityExt(lambda, m, pp) * ds;
            }
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauinteract < tau)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauinteract, tau0, tau, s0, s));
                return true;
            }
        }
    }

    // the interaction point is outside of the path
    return false;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::setInteractionPointUsingScatteringAndAbsorption(PhotonPacket* pp, double tauinteract) const
{
    auto generator = getPathSegmentGenerator(_grid, pp);
    double tauSca = 0.;
    double tauAbs = 0.;
    double s = 0.;

    // loop over the segments of the path until the interaction optical depth is reached or the path ends

    // --> single medium, spatially constant cross sections
    if (_config->hasSingleConstantSectionMedium())
    {
        double sectionSca = mix(0, 0)->sectionSca(pp->wavelength());
        double sectionAbs = mix(0, 0)->sectionAbs(pp->wavelength());
        while (generator->next())
        {
            // remember the cumulative scattering optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tauSca0 = tauSca;
            double s0 = s;

            // calculate the cumulative optical depths and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0) tauSca += sectionSca * _state.numberDensity(m, 0) * generator->ds();
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauinteract < tauSca)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauinteract, tauSca0, tauSca, s0, s),
                                        tauinteract * sectionAbs / sectionSca);
                return true;
            }
        }
    }

    // --> multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantSectionMedia())
    {
        ShortArray sectionScav(_numMedia);
        ShortArray sectionAbsv(_numMedia);
        for (int h = 0; h != _numMedia; ++h)
        {
            sectionScav[h] = mix(0, h)->sectionSca(pp->wavelength());
            sectionAbsv[h] = mix(0, h)->sectionAbs(pp->wavelength());
        }
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tauSca0 = tauSca;
            double tauAbs0 = tauAbs;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                for (int h = 0; h != _numMedia; ++h)
                {
                    double ns = _state.numberDensity(m, h) * ds;
                    tauSca += sectionScav[h] * ns;
                    tauAbs += sectionAbsv[h] * ns;
                }
            }
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauinteract < tauSca)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauinteract, tauSca0, tauSca, s0, s),
                                        NR::interpolateLinLin(tauinteract, tauSca0, tauSca, tauAbs0, tauAbs));
                return true;
            }
        }
    }

    // --> spatially variable cross sections
    else
    {
        while (generator->next())
        {
            // remember the cumulative optical depth and distance at the start of this segment
            // so that we can interpolate the interaction point should it happen to be inside this segment
            double tauSca0 = tauSca;
            double tauAbs0 = tauAbs;
            double s0 = s;

            // calculate the cumulative optical depth and distance at the end of this segment
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                double lambda = pp->perceivedWavelength(_state.bulkVelocity(m), _config->hubbleExpansionRate() * s);
                tauSca += opacitySca(lambda, m, pp) * ds;
                tauAbs += opacityAbs(lambda, m, pp) * ds;
            }
            s += ds;

            // if the interaction point is inside this segment, store it in the photon packet
            if (tauinteract < tauSca)
            {
                pp->setInteractionPoint(m, NR::interpolateLinLin(tauinteract, tauSca0, tauSca, s0, s),
                                        NR::interpolateLinLin(tauinteract, tauSca0, tauSca, tauAbs0, tauAbs));
                return true;
            }
        }
    }

    // the interaction point is outside of the path
    return false;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::getExtinctionOpticalDepth(const PhotonPacket* pp, double distance) const
{
    // abort if the packet's contribution is zero to begin with
    double L = pp->luminosity();
    if (L <= 0) return std::numeric_limits<double>::infinity();

    // if extinction is always positive, determine the optical depth at which the packet's contribution becomes zero
    double taumax = _config->hasNegativeExtinction() ? std::numeric_limits<double>::infinity() : std::log(L) + 745;

    // determine the geometric details of the path and calculate the optical depth at the same time
    auto generator = getPathSegmentGenerator(_grid, pp);
    double tau = 0.;
    double s = 0.;

    // single medium, spatially constant cross sections
    if (_config->hasSingleConstantSectionMedium())
    {
        double section = mix(0, 0)->sectionExt(pp->wavelength());
        while (generator->next())
        {
            if (generator->m() >= 0)
            {
                tau += section * _state.numberDensity(generator->m(), 0) * generator->ds();
                if (tau >= taumax) return std::numeric_limits<double>::infinity();
            }
            s += generator->ds();
            if (s > distance) break;
        }
    }

    // multiple media, spatially constant cross sections
    else if (_config->hasMultipleConstantSectionMedia())
    {
        ShortArray sectionv(_numMedia);
        for (int h = 0; h != _numMedia; ++h) sectionv[h] = mix(0, h)->sectionExt(pp->wavelength());
        while (generator->next())
        {
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                for (int h = 0; h != _numMedia; ++h) tau += sectionv[h] * _state.numberDensity(m, h) * ds;
                if (tau >= taumax) return std::numeric_limits<double>::infinity();
            }
            s += ds;
            if (s > distance) break;
        }
    }

    // spatially variable cross sections
    else
    {
        while (generator->next())
        {
            double ds = generator->ds();
            int m = generator->m();
            if (m >= 0)
            {
                double lambda = pp->perceivedWavelength(_state.bulkVelocity(m), _config->hubbleExpansionRate() * s);
                tau += opacityExt(lambda, m, pp) * ds;
                if (tau >= taumax) return std::numeric_limits<double>::infinity();
            }
            s += ds;
            if (s > distance) break;
        }
    }

    return tau;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::getExtinctionOpticalDepth(const SpatialGridPath* path, double lambda,
                                               MaterialMix::MaterialType type) const
{
    // determine the geometric details of the path and calculate the optical depth at the same time
    auto generator = getPathSegmentGenerator(_grid, path);
    double tau = 0.;
    while (generator->next())
    {
        if (generator->m() >= 0) tau += opacityExt(lambda, generator->m(), type) * generator->ds();
    }
    return tau;
}

////////////////////////////////////////////////////////////////////

void MediumSystem::clearRadiationField(bool primary)
{
    if (primary)
    {
        _rf1.setToZero();
        std::fill(_labs1.begin(), _labs1.end(), 0.);
        _labsOutside1 = 0.;
    }
    else
    {
        _rf2c.setToZero();
        std::fill(_labs2c.begin(), _labs2c.end(), 0.);
        _labsOutside2c = 0.;
    }
}

////////////////////////////////////////////////////////////////////

// This helper class calculates the absorption opacity of the dust media in a given cell, at the specified
// wavelength, for the specified photon packet traveling along its current path. If the dust absorption cross
// sections are spatially constant, the constructor obtains these cross sections for each dust medium component,
// so that the work that does not depend on the cell is performed once per path rather than once per segment.
// Otherwise, the opacity is obtained from the material mix of each dust component in the cell, which handles
// material mix families, extra specific state variables, and polarization-dependent absorption.
class MediumSystem::DustAbsorptionOpacity
{
public:
    DustAbsorptionOpacity(const MediumSystem* ms, double lambda, const PhotonPacket* pp)
        : _ms{ms}, _lambda{lambda}, _pp{pp}, _sectionv(ms->_dust_hv.size())
    {
        if (_ms->hasExactDustAbsorption() && _ms->_config->hasConstantDustAbsorptionSections())
        {
            int numDust = _ms->_dust_hv.size();
            for (int i = 0; i != numDust; ++i) _sectionv[i] = _ms->mix(0, _ms->_dust_hv[i])->sectionAbs(_lambda);
            _constant = true;
        }
    }

    // returns false if the simulation has no dust emission, so that no opacity is needed
    bool active() const { return _ms->hasExactDustAbsorption(); }

    // returns the absorption opacity of the dust media in cell m
    double opacity(int m) const
    {
        double result = 0.;
        int numDust = _ms->_dust_hv.size();
        if (_constant)
        {
            for (int i = 0; i != numDust; ++i) result += _sectionv[i] * _ms->_state.numberDensity(m, _ms->_dust_hv[i]);
        }
        else
        {
            for (int h : _ms->_dust_hv)
            {
                MaterialState mst(_ms->_state, m, h);
                result += _ms->mix(m, h)->opacityAbs(_lambda, &mst, _pp);
            }
        }
        return result;
    }

private:
    const MediumSystem* _ms;
    double _lambda;
    const PhotonPacket* _pp;
    bool _constant{false};
    ShortArray _sectionv;  // the absorption cross section of each dust component, if constant
};

////////////////////////////////////////////////////////////////////

void MediumSystem::storeRadiationField(bool primary, const PhotonPacket* pp)
{
    // the targets for the radiation field and the dust-absorbed luminosity
    CellTable& rf = primary ? _rf1 : _rf2c;
    vector<double>& labs = primary ? _labs1 : _labs2c;
    double LabsOutside = 0.;  // the part of the absorbed luminosity along the path outside the wavelength grid

    // the rebalance tallies for this path, which are active only for dust emission during a rebalanced iteration
    // (in which case the dust absorption opacity is always active as well); the launched luminosity is added for the
    // first path of the photon packet only
    SelfHeatingRebalancer::Tally tally;
    if (_rebalancer && !primary) tally = _rebalancer->tally(pp->emissionCell());
    if (tally && pp->numScatt() == 0) tally.launch(pp->luminosity());

    // use a faster version in case there are no kinematics
    if (_config->hasConstantPerceivedWavelength())
    {
        int ell = _wavelengthGrid->bin(pp->wavelength());
        DustAbsorptionOpacity dust(this, pp->wavelength(), pp);
        if (ell < 0 && !dust.active()) return;

        double luminosity = pp->luminosity();
        double lnExtBeg = 0.;  // extinction factor and its logarithm at begin of current segment
        double extBeg = 1.;
        for (const auto& segment : pp->segments())
        {
            double lnExtEnd = -segment.tauExt();  // extinction factor and its logarithm at end of current segment
            double extEnd = exp(lnExtEnd);
            int m = segment.m();
            if (m >= 0)
            {
                // use this flavor of the lnmean function to avoid recalculating the logarithm of the extinction
                double extMean = SpecialFunctions::lnmean(extEnd, extBeg, lnExtEnd, lnExtBeg);
                double Lds = luminosity * extMean * segment.ds();
                if (ell >= 0) LockFree::add(rf(m, ell), Lds);
                if (dust.active())
                {
                    double Labs = dust.opacity(m) * Lds;
                    LockFree::add(labs[m], Labs);
                    if (ell < 0) LabsOutside += Labs;
                    if (tally) tally.absorb(m, Labs);
                }
            }
            lnExtBeg = lnExtEnd;
            extBeg = extEnd;
        }
        if (tally) tally.escape(luminosity * extBeg);
    }
    else
    {
        double lnExtBeg = 0.;  // extinction factor and its logarithm at begin of current segment
        double extBeg = 1.;
        for (const auto& segment : pp->segments())
        {
            double lnExtEnd = -segment.tauExt();  // extinction factor and its logarithm at end of current segment
            double extEnd = exp(lnExtEnd);
            int m = segment.m();
            if (m >= 0)
            {
                double lambda =
                    pp->perceivedWavelength(_state.bulkVelocity(m), _config->hubbleExpansionRate() * segment.s());
                int ell = _wavelengthGrid->bin(lambda);
                DustAbsorptionOpacity dust(this, lambda, pp);
                if (ell >= 0 || dust.active())
                {
                    // use this flavor of the lnmean function to avoid recalculating the logarithm of the extinction
                    double extMean = SpecialFunctions::lnmean(extEnd, extBeg, lnExtEnd, lnExtBeg);
                    double Lds = pp->perceivedLuminosity(lambda) * extMean * segment.ds();
                    if (ell >= 0) LockFree::add(rf(m, ell), Lds);
                    if (dust.active())
                    {
                        double Labs = dust.opacity(m) * Lds;
                        LockFree::add(labs[m], Labs);
                        if (ell < 0) LabsOutside += Labs;
                        if (tally) tally.absorb(m, Labs);
                    }
                }
            }
            lnExtBeg = lnExtEnd;
            extBeg = extEnd;
        }
        if (tally) tally.escape(pp->luminosity() * extBeg);
    }

    // add the part absorbed outside the grid once per path, so that threads do not contend for it in the loop
    if (LabsOutside > 0.) LockFree::add(primary ? _labsOutside1 : _labsOutside2c, LabsOutside);
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the specified part of the sum of the specified array, after summing that part over all processes,
    // or zero if the sum is zero
    double fractionOfSum(double part, const vector<double>& array)
    {
        Array partv(part, 1);
        ProcessManager::sumToAll(partv);
        double sum = std::accumulate(array.begin(), array.end(), 0.);
        return sum > 0. ? partv[0] / sum : 0.;
    }
}

////////////////////////////////////////////////////////////////////

void MediumSystem::communicateRadiationField(bool primary)
{
    if (primary)
    {
        ProcessManager::sumToAll(_rf1.data());
        if (_labs1.size())
        {
            ProcessManager::sumToAll(_labs1);
            _fractionOutside1 = fractionOfSum(_labsOutside1, _labs1);
        }
    }
    else
    {
        ProcessManager::sumToAll(_rf2c.data());
        _rf2 = _rf2c;
        if (_labs2c.size())
        {
            ProcessManager::sumToAll(_labs2c);
            _labs2 = _labs2c;
            _fractionOutside2 = fractionOfSum(_labsOutside2c, _labs2);
        }
    }
}

////////////////////////////////////////////////////////////////////

std::pair<double, double> MediumSystem::totalDustAbsorbedLuminosity() const
{
    // the arrays are identical on all processes after communicateRadiationField(), so no communication is needed
    double Labs1 = std::accumulate(_labs1.begin(), _labs1.end(), 0.);
    double Labs2 = std::accumulate(_labs2.begin(), _labs2.end(), 0.);
    return std::make_pair(Labs1, Labs2);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::dustAbsorbedFractionOutsideGrid(bool primary) const
{
    return primary ? _fractionOutside1 : _fractionOutside2;
}

////////////////////////////////////////////////////////////////////

void MediumSystem::beginRebalance()
{
    // calculate the dust mass in each cell
    vector<double> dustMassv(_numCells);
    for (int m = 0; m != _numCells; ++m)
    {
        double rho = 0.;
        for (int h : _dust_hv) rho += massDensity(m, h);
        dustMassv[m] = rho * volume(m);
    }

    _rebalancer->begin(_labs1, _labs2, dustMassv);
}

////////////////////////////////////////////////////////////////////

void MediumSystem::endRebalance()
{
    _rebalancer->end(_labs1, _labs2);

    // rescale the stable secondary radiation field and absorbed luminosity; the factors are identical on all processes
    int numBins = _wavelengthGrid->numBins();
    for (int m = 0; m != _numCells; ++m)
    {
        double f = _rebalancer->factor(m);
        if (f != 1.)
        {
            _labs2[m] *= f;
            for (int ell = 0; ell != numBins; ++ell) _rf2(m, ell) *= f;
        }
    }
}

////////////////////////////////////////////////////////////////////

double MediumSystem::radiationField(int m, int ell) const
{
    double rf = 0.;
    if (_rf1.size()) rf += _rf1(m, ell);
    if (_rf2.size()) rf += _rf2(m, ell);
    return rf;
}

////////////////////////////////////////////////////////////////////

Array MediumSystem::meanIntensity(int m) const
{
    int numWavelengths = _wavelengthGrid->numBins();
    Array Jv(numWavelengths);
    double factor = 1. / (4. * M_PI * _state.volume(m));
    for (int ell = 0; ell < numWavelengths; ell++)
    {
        Jv[ell] = radiationField(m, ell) * factor / _wavelengthGrid->effectiveWidth(ell);
    }
    return Jv;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeTemperature(int m, int h) const
{
    // get the radiation field, if available
    Array Jv;
    if (_config->hasRadiationField()) Jv = meanIntensity(m);

    // obtain the temperature for the requested component and cell
    if (massDensity(m, h) > 0.)
    {
        MaterialState mst(_state, m, h);
        return mix(m, h)->indicativeTemperature(&mst, Jv);
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeTemperature(int m, MaterialMix::MaterialType type) const
{
    // get the radiation field, if available
    Array Jv;
    if (_config->hasRadiationField()) Jv = meanIntensity(m);

    // obtain the temperature and weight for each component of the requested type
    double sumRhoT = 0.;
    double sumRho = 0.;
    for (int h = 0; h != _numMedia; ++h)
    {
        if (isMaterialType(type, h))
        {
            double rho = massDensity(m, h);
            if (rho > 0.)
            {
                MaterialState mst(_state, m, h);
                double T = mix(m, h)->indicativeTemperature(&mst, Jv);
                sumRhoT += rho * T;
                sumRho += rho;
            }
        }
    }

    // return the average, if there is one
    return sumRho > 0. ? sumRhoT / sumRho : 0.;
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeDustTemperature(int m) const
{
    return indicativeTemperature(m, MaterialMix::MaterialType::Dust);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeElectronTemperature(int m) const
{
    return indicativeTemperature(m, MaterialMix::MaterialType::Electrons);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::indicativeGasTemperature(int m) const
{
    return indicativeTemperature(m, MaterialMix::MaterialType::Gas);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::dustLuminosity(int m) const
{
    double Labs = 0.;
    if (_labs1.size()) Labs += _labs1[m];
    if (_labs2.size()) Labs += _labs2[m];
    return Labs;
}

////////////////////////////////////////////////////////////////////

Array MediumSystem::dustEmissionSpectrum(int m) const
{
    const Array& Jv = meanIntensity(m);
    Array ev(_numDustEmissionWavelengths);
    for (int h : _dust_hv)
    {
        MaterialState mst(_state, m, h);
        ev += mix(m, h)->emissionSpectrum(&mst, Jv);
    }
    return ev;
}

////////////////////////////////////////////////////////////////////

Array MediumSystem::continuumEmissionSpectrum(int m, int h) const
{
    const Array& Jv = meanIntensity(m);
    MaterialState mst(_state, m, h);
    return mix(m, h)->emissionSpectrum(&mst, Jv);
}

////////////////////////////////////////////////////////////////////

Array MediumSystem::lineEmissionSpectrum(int m, int h) const
{
    const Array& Jv = meanIntensity(m);
    MaterialState mst(_state, m, h);
    return mix(m, h)->lineEmissionSpectrum(&mst, Jv);
}

////////////////////////////////////////////////////////////////////

void MediumSystem::recordAggregates()
{
    // collect the aggregate series and the (medium component index, custom index) pair for each of them
    auto aggregates = _history->aggregates();
    if (aggregates.empty()) return;
    vector<std::pair<int, int>> variables;
    variables.reserve(aggregates.size());
    for (auto series : aggregates)
    {
        auto rule = series->aggregationRule();
        auto position = std::find(_media.begin(), _media.end(), rule->medium);
        if (position == _media.end()) throw FATALERROR("An aggregate series refers to an unknown medium component");
        variables.emplace_back(static_cast<int>(position - _media.begin()), rule->customIndex);
    }

    // calculate the volume integrals in a single pass over the cells, and record them
    auto integrals = _state.volumeIntegrals(variables);
    for (size_t k = 0; k != aggregates.size(); ++k) aggregates[k]->set(integrals[k]);
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::updateDynamicStateRecipes()
{
    auto log = find<Log>();
    auto parfac = find<ParallelFactory>();
    auto& recipes = dynamicStateOptions()->recipes();

    // tell all recipes to begin the update cycle
    for (auto recipe : recipes) recipe->beginUpdate(_numCells);

    // update status for each cell
    std::vector<UpdateStatus> flags(_numCells);

    // loop over the spatial cells in parallel
    log->info("Updating medium state for " + std::to_string(_numCells) + " cells...");
    log->infoSetElapsed(_numCells);
    parfac->parallelDistributed()->call(_numCells, [this, log, &recipes, &flags](size_t firstIndex, size_t numIndices) {
        while (numIndices)
        {
            size_t currentChunkSize = min(logProgressChunkSize, numIndices);
            for (size_t m = firstIndex; m != firstIndex + currentChunkSize; ++m)
            {
                const Array& Jv = meanIntensity(m);

                // tell all recipes to update the state for all media in this cell
                for (auto recipe : recipes)
                {
                    for (int h = 0; h != _numMedia; ++h)
                    {
                        MaterialState mst(_state, m, h);
                        flags[m].update(recipe->update(&mst, Jv));
                    }
                }
            }
            log->infoIfElapsed("Updated medium state: ", currentChunkSize);
            firstIndex += currentChunkSize;
            numIndices -= currentChunkSize;
        }
    });

    // synchronize the updated state between processes
    auto [numUpdated, numNotConverged] = _state.synchronize(flags);

    // log statistics
    log->info("  Updated cells: " + std::to_string(numUpdated) + " out of " + std::to_string(_numCells) + " ("
              + StringUtils::toString(100. * numUpdated / _numCells, 'f', 2) + " %)");
    log->info("  Not converged: " + std::to_string(numNotConverged) + " out of " + std::to_string(_numCells) + " ("
              + StringUtils::toString(100. * numNotConverged / _numCells, 'f', 2) + " %)");

    // record the aggregates of the updated state; a recipe that needs aggregate information declares its own
    // aggregate series in the iteration history during setup, and retrieves it in endUpdate()
    recordAggregates();

    // tell all recipes to end the update cycle and collect convergence info
    bool converged = true;
    for (auto recipe : recipes) converged &= recipe->endUpdate(_numCells, numUpdated, numNotConverged);
    return converged;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::updateDynamicStateMedia(bool primary)
{
    auto log = find<Log>();
    auto parfac = find<ParallelFactory>();
    const vector<int>& hv = primary ? _pdms_hv : _sdms_hv;

    // determine the medium components that update their state in stages
    vector<int> stagedHv;
    for (int h : hv)
        if (_media[h]->capabilities().stagedStateUpdate) stagedHv.push_back(h);

    // for these components, prepare the update of each cell (stage 1) in parallel, synchronize the prepared state
    // between processes, and begin the update (stage 2) on each process
    if (!stagedHv.empty())
    {
        std::vector<UpdateStatus> flags(_numCells);
        log->info("Preparing medium state update for " + std::to_string(_numCells) + " cells...");
        log->infoSetElapsed(_numCells);
        parfac->parallelDistributed()->call(
            _numCells, [this, log, &stagedHv, &flags](size_t firstIndex, size_t numIndices) {
                while (numIndices)
                {
                    size_t currentChunkSize = min(logProgressChunkSize, numIndices);
                    for (size_t m = firstIndex; m != firstIndex + currentChunkSize; ++m)
                    {
                        const Array& Jv = meanIntensity(m);
                        for (int h : stagedHv)
                        {
                            MaterialState mst(_state, m, h);
                            flags[m].update(mix(m, h)->prepareSpecificStateUpdate(&mst, Jv));
                        }
                    }
                    log->infoIfElapsed("Prepared medium state update: ", currentChunkSize);
                    firstIndex += currentChunkSize;
                    numIndices -= currentChunkSize;
                }
            });
        _state.synchronize(flags);
        for (int h : stagedHv) mix(0, h)->beginSpecificStateUpdate(_state, h);
    }

    // update status for each cell
    std::vector<UpdateStatus> flags(_numCells);

    // loop over the spatial cells in parallel (stage 3 for the components that update their state in stages)
    log->info("Updating medium state for " + std::to_string(_numCells) + " cells...");
    log->infoSetElapsed(_numCells);
    parfac->parallelDistributed()->call(_numCells, [this, &hv, log, &flags](size_t firstIndex, size_t numIndices) {
        while (numIndices)
        {
            size_t currentChunkSize = min(logProgressChunkSize, numIndices);
            for (size_t m = firstIndex; m != firstIndex + currentChunkSize; ++m)
            {
                const Array& Jv = meanIntensity(m);
                for (int h : hv)
                {
                    MaterialState mst(_state, m, h);
                    flags[m].update(mix(m, h)->updateSpecificState(&mst, Jv));
                }
            }
            log->infoIfElapsed("Updated medium state: ", currentChunkSize);
            firstIndex += currentChunkSize;
            numIndices -= currentChunkSize;
        }
    });

    // synchronize the updated state between processes
    auto [numUpdated, numNotConverged] = _state.synchronize(flags);

    // log statistics
    log->info("  Updated cells: " + std::to_string(numUpdated) + " out of " + std::to_string(_numCells) + " ("
              + StringUtils::toString(100. * numUpdated / _numCells, 'f', 2) + " %)");
    if (numNotConverged)
        log->info("  Not converged: " + std::to_string(numNotConverged) + " out of " + std::to_string(_numCells) + " ("
                  + StringUtils::toString(100. * numNotConverged / _numCells, 'f', 2) + " %)");

    // end the update (stage 4) on each process for the components that update their state in stages
    for (int h : stagedHv) mix(0, h)->endSpecificStateUpdate(_state, h);

    // record the aggregates of the updated state
    recordAggregates();

    // collect convergence info from the material mixes, passing each a history scope bound to its medium component
    bool converged = true;
    for (int h : hv)
    {
        HistoryScope scope(_history, _media[h]);
        converged &= mix(0, h)->isSpecificStateConverged(_numCells, numUpdated, numNotConverged, scope);
    }
    return converged;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::updatePrimaryDynamicMediumState()
{
    bool converged = true;
    if (_config->hasDynamicStateRecipes()) converged &= updateDynamicStateRecipes();
    if (_config->hasPrimaryDynamicStateMedia()) converged &= updateDynamicStateMedia(true);
    return converged;
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::updateSecondaryDynamicMediumState()
{
    bool converged = true;
    if (_config->hasSecondaryDynamicStateMedia()) converged &= updateDynamicStateMedia(false);
    return converged;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // the ids of the iteration history series and cell windows kept by the medium system for dynamic refinement,
    // keyed on the dynamic refinement options, on a refinement criterion, or on a refinement condition
    enum RefinementId {
        Subdivided,
        Rounds,
        Accumulated,
        NumCells,
        PacketsPerCell,
        CellMemory,
        Settled,
        Candidates,
        LowestPriority,
        CriterionCandidates,
        BlockedByLevel,
        BlockedByColumnDensity,
        Field
    };

    // the bit flags recording the outcome of the evaluation of a cell for a refinement criterion
    constexpr int candidateFlag = 1;      // the cell is a candidate for the criterion
    constexpr int levelFlag = 2;          // the criterion flags the cell, but it has reached the level limit
    constexpr int columnDensityFlag = 4;  // the criterion flags the cell, but it has reached the column density limit

    // the number of photon packets per cell launched in an iteration below which a warning is issued
    constexpr double minPacketsPerCell = 100.;

    // returns a description of the specified refinement criterion for use in log messages and series descriptions
    string criterionDescription(const RefinementCriterion* criterion)
    {
        return criterion->type() + " on " + criterion->fieldDescription();
    }

    // returns a description of the specified iteration loop for use in log messages
    string loopDescription(IterationHistory::Loop loop)
    {
        switch (loop)
        {
            case IterationHistory::Loop::None: break;
            case IterationHistory::Loop::Primary: return "primary emission";
            case IterationHistory::Loop::Secondary: return "secondary emission";
            case IterationHistory::Loop::Merged: return "merged primary and secondary emission";
        }
        return "unknown";
    }
}

////////////////////////////////////////////////////////////////////

void MediumSystem::setupDynamicRefinement()
{
    // dynamic refinement requires a tree grid with refinement options in a simulation that iterates
    _tree = dynamic_cast<TreeSpatialGrid*>(_grid);
    if (!_tree || !_tree->dynamicRefinementOptions()
        || (!_config->hasPrimaryIterations() && !_config->hasSecondaryIterations()))
    {
        _tree = nullptr;
        return;
    }
    _refinement = _tree->dynamicRefinementOptions();
    auto log = find<Log>();

    // determine the iteration loops of the simulation and the number of photon packets launched in each iteration
    using Loop = IterationHistory::Loop;
    vector<std::pair<Loop, double>> loops;
    double numPrimary = _config->numPrimaryIterationPackets();
    double numSecondary = _config->numSecondaryIterationPackets();
    if (_config->hasPrimaryIterations()) loops.emplace_back(Loop::Primary, numPrimary);
    if (_config->hasMergedIterations())
        loops.emplace_back(Loop::Merged, numPrimary + numSecondary);
    else if (_config->hasSecondaryIterations())
        loops.emplace_back(Loop::Secondary, numSecondary);

    // keep the loops in which at least one criterion takes part
    vector<std::pair<Loop, double>> refiningLoops;
    for (const auto& loop : loops)
    {
        for (auto criterion : _refinement->criteria())
        {
            if (criterion->changesIn(loop.first))
            {
                refiningLoops.push_back(loop);
                break;
            }
        }
    }
    if (refiningLoops.empty())
        log->warning("None of the dynamic grid refinement criteria takes part in an iteration loop of this simulation");

    // the remainder applies only if there is a cap on the number of cells
    int maxCells = _refinement->maxCells();
    if (!maxCells) return;
    if (maxCells <= _numCells)
    {
        log->warning("The maximum number of cells for dynamic grid refinement (" + std::to_string(maxCells)
                     + ") does not exceed the number of cells in the initial grid (" + std::to_string(_numCells) + ")");
        return;
    }

    // reserve capacity for the maximum number of cells in all per-cell data structures
    _tree->reserveCells(maxCells);
    _state.reserveCells(maxCells);
    for (auto table : {&_rf1, &_rf2, &_rf2c})
        if (table->size()) table->reserveCells(maxCells);
    for (auto labs : {&_labs1, &_labs2, &_labs2c})
        if (labs->size()) labs->reserve(maxCells);
    if (_mixPerCell) _mixv.reserve(static_cast<size_t>(maxCells) * _numMedia);
    _history->reserveCells(maxCells);

    // log the reserved memory, including the cell windows of the refinement criteria and conditions declared later
    size_t numFields = 0;
    for (auto criterion : _refinement->criteria()) numFields += 1 + criterion->conditions().size();
    double perCell = bytesPerCell() + numFields * sizeof(double);
    log->info("Reserved memory for the per-cell data of up to " + std::to_string(maxCells)
              + " cells: " + StringUtils::toMemSizeString(perCell * maxCells) + " per process ("
              + StringUtils::toString(perCell, 'f', 0) + " bytes per cell)");

    // warn if the memory of the process with all of that capacity in use would exceed the physical memory
    double needed = System::currentMemoryUsage() + perCell * (maxCells - _numCells);
    size_t available = System::availableMemory();
    if (available && needed > available)
        log->warning("The memory needed per process with the maximum number of cells ("
                     + StringUtils::toMemSizeString(needed) + ") exceeds the physical memory of the host ("
                     + StringUtils::toMemSizeString(available) + ")");

    // warn if the number of photon packets per cell would be low for a grid with the maximum number of cells
    for (const auto& loop : refiningLoops)
    {
        double packetsPerCell = loop.second / maxCells;
        if (packetsPerCell < minPacketsPerCell)
            log->warning("With the maximum number of cells, each " + loopDescription(loop.first)
                         + " iteration launches only " + StringUtils::toString(packetsPerCell, 'g', 3)
                         + " photon packets per cell");
    }
}

////////////////////////////////////////////////////////////////////

double MediumSystem::bytesPerCell() const
{
    double bytes = _state.bytesPerCell();
    bytes += _rf1.bytesPerCell() + _rf2.bytesPerCell() + _rf2c.bytesPerCell();
    for (auto labs : {&_labs1, &_labs2, &_labs2c})
        if (labs->size()) bytes += sizeof(double);
    if (_mixPerCell) bytes += _numMedia * sizeof(const MaterialMix*);
    bytes += _history->allCellWindows().size() * sizeof(double);
    if (_tree) bytes += _tree->bytesPerAddedCell();
    return bytes;
}

////////////////////////////////////////////////////////////////////

void MediumSystem::growCellData(const vector<int>& subdividedCells, const vector<int>& parentCells)
{
    int numOldCells = _numCells;
    _numCells = _grid->numCells();

    // copy the values of the parent to each new cell
    _state.appendCells(parentCells);
    for (auto table : {&_rf1, &_rf2, &_rf2c})
        if (table->size()) table->appendCells(parentCells);
    for (auto labs : {&_labs1, &_labs2, &_labs2c})
        if (labs->size())
            for (int parent : parentCells) labs->push_back((*labs)[parent]);
    if (_mixPerCell)
        for (int parent : parentCells)
            for (int h = 0; h != _numMedia; ++h) _mixv.push_back(_mixv[parent * _numMedia + h]);
    _history->appendCells(parentCells);

    // set the volume of each child, which still holds the volume of its parent, and divide the radiation field, the
    // dust-absorbed luminosities, and the extensive state variables in proportion to the volume of the children
    auto resize = [this](int m) {
        double volume = _grid->volume(m);
        double oldVolume = _state.volume(m);
        double fraction = oldVolume > 0. ? volume / oldVolume : 0.;
        for (auto table : {&_rf1, &_rf2, &_rf2c})
            if (table->size()) table->scaleCell(m, fraction);
        for (auto labs : {&_labs1, &_labs2, &_labs2c})
            if (labs->size()) (*labs)[m] *= fraction;
        _state.resizeCell(m, volume);
    };
    for (int m : subdividedCells) resize(m);
    for (int m = numOldCells; m != _numCells; ++m) resize(m);

    // verify that each per-cell data structure has grown to the new number of cells
    auto verify = [this](size_t size, size_t expected, string name) {
        if (size != expected)
            throw FATALERROR("The size of the " + name + " (" + std::to_string(size)
                             + ") does not match the number of cells after dynamic refinement ("
                             + std::to_string(_numCells) + ")");
    };
    verify(_state.numCells(), _numCells, "medium state");
    for (auto table : {&_rf1, &_rf2, &_rf2c})
        if (table->size()) verify(table->numCells(), _numCells, "radiation field table");
    for (auto labs : {&_labs1, &_labs2, &_labs2c})
        if (labs->size()) verify(labs->size(), _numCells, "dust-absorbed luminosity array");
    if (_mixPerCell) verify(_mixv.size(), static_cast<size_t>(_numCells) * _numMedia, "material mix array");
    for (auto window : _history->allCellWindows()) verify(window->numCells(), _numCells, "cell window");
}

////////////////////////////////////////////////////////////////////

bool MediumSystem::updateDynamicRefinement(double numPackets)
{
    if (!_refinement) return true;

    // determine the criteria that take part in the current iteration loop
    auto loop = _history->loop();
    vector<const RefinementCriterion*> criteria;
    for (auto criterion : _refinement->criteria())
        if (criterion->changesIn(loop)) criteria.push_back(criterion);
    if (criteria.empty()) return true;
    int numCriteria = criteria.size();

    auto log = find<Log>();
    auto options = _refinement;
    auto lifetime = IterationHistory::Lifetime::Loop;
    int numLeftOut = options->numIterationsAfterRefinement();
    int numAveraged = options->numAveragedIterations();
    int iteration = _history->iteration();
    int maxIterations =
        loop == IterationHistory::Loop::Primary ? _config->maxPrimaryIterations() : _config->maxSecondaryIterations();

    // declare all scalar series in every iteration (declaring is idempotent), so that they exist from the first
    // iteration of the loop, including those that receive a value only in decision iterations
    auto& subdivided =
        _history->scalarSeries(options, Subdivided, numLeftOut + 1, lifetime, "number of subdivided cells");
    auto& rounds = _history->scalarSeries(options, Rounds, 2, lifetime, "number of refinement rounds");
    auto& accumulated =
        _history->scalarSeries(options, Accumulated, 1, lifetime, "number of iterations in the refinement window");
    auto& cells = _history->scalarSeries(options, NumCells, 1, lifetime, "total number of cells");
    auto& packets = _history->scalarSeries(options, PacketsPerCell, 1, lifetime, "number of photon packets per cell");
    auto& memory = _history->scalarSeries(options, CellMemory, 1, lifetime, "memory held by the per-cell data (bytes)");
    auto& settled = _history->scalarSeries(options, Settled, 1, lifetime, "refinement settled");
    auto& candidates = _history->scalarSeries(options, Candidates, 1, lifetime, "number of refinement candidates");
    auto& lowest =
        _history->scalarSeries(options, LowestPriority, 1, lifetime, "lowest priority among the subdivided cells");
    vector<ScalarSeries*> criterionCandidates, blockedByLevel, blockedByColumnDensity;
    for (auto criterion : criteria)
    {
        string description = criterionDescription(criterion);
        criterionCandidates.push_back(&_history->scalarSeries(criterion, CriterionCandidates, 1, lifetime,
                                                              "number of candidates for " + description));
        blockedByLevel.push_back(&_history->scalarSeries(
            criterion, BlockedByLevel, 1, lifetime, "number of cells blocked by the level limit of " + description));
        blockedByColumnDensity.push_back(
            &_history->scalarSeries(criterion, BlockedByColumnDensity, 1, lifetime,
                                    "number of cells blocked by the column density limit of " + description));
    }

    // record the values set in every iteration; the number of rounds is carried over from the previous iteration
    subdivided.set(0.);
    int numRounds = rounds.has(1) ? static_cast<int>(rounds.value(1)) : 0;
    rounds.set(numRounds);
    cells.set(_numCells);
    packets.set(numPackets / _numCells);
    memory.set(bytesPerCell() * _numCells);

    // this function records whether the refinement has settled, warns if the loop ends before it has settled, and
    // returns the result
    auto finish = [&settled, &log, iteration, maxIterations, loop](bool hasSettled) {
        settled.set(hasSettled ? 1. : 0.);
        if (!hasSettled && iteration >= maxIterations)
            log->warning("The " + loopDescription(loop)
                         + " iteration loop ends before the dynamic grid refinement has settled");
        return hasSettled;
    };

    // during the initial iterations, refinement is not considered at all
    if (iteration <= options->numIterationsBeforeRefinement())
    {
        accumulated.set(0.);
        return finish(false);
    }

    // an iteration shortly after a refinement round is left out of the averaging
    bool leftOut = false;
    for (int lag = 1; lag <= numLeftOut; ++lag) leftOut |= subdivided.has(lag) && subdivided.value(lag) > 0.;

    // obtain the cell windows of the criteria and of their conditions, which share a single schedule
    vector<CellWindow*> fieldWindows;
    vector<vector<std::pair<const RefinementCondition*, CellWindow*>>> conditionWindows(numCriteria);
    for (int c = 0; c != numCriteria; ++c)
    {
        auto criterion = criteria[c];
        fieldWindows.push_back(&_history->cellWindow(criterion, Field, _numCells, lifetime,
                                                     "window-averaged " + criterion->fieldDescription()));
        for (auto condition : criterion->conditions())
            conditionWindows[c].emplace_back(condition,
                                             &_history->cellWindow(condition, Field, _numCells, lifetime,
                                                                   "window-averaged " + condition->description()));
    }

    // accumulate the current values, unless the iteration is left out
    if (leftOut)
    {
        log->info("Dynamic grid refinement: iteration left out of the averaging after a refinement round");
    }
    else
    {
        auto parallel = find<ParallelFactory>()->parallelDistributed();
        for (int c = 0; c != numCriteria; ++c)
        {
            auto criterion = criteria[c];
            if (criterion->hasCostlyValue())
            {
                // calculate the values in parallel, each cell on a single process, and sum them over the processes
                Array values(_numCells);
                parallel->call(_numCells, [criterion, &values](size_t firstIndex, size_t numIndices) {
                    for (size_t m = firstIndex; m != firstIndex + numIndices; ++m) values[m] = criterion->value(m);
                });
                ProcessManager::sumToAll(values);
                fieldWindows[c]->accumulate([&values](int m) { return values[m]; });
            }
            else
            {
                fieldWindows[c]->accumulate([criterion](int m) { return criterion->value(m); });
            }
            for (auto& entry : conditionWindows[c])
            {
                auto condition = entry.first;
                entry.second->accumulate([condition](int m) { return condition->value(m); });
            }
        }
    }
    int numAccumulated = fieldWindows[0]->numIterations();
    accumulated.set(numAccumulated);

    // if the windows are not yet complete, the refinement has not settled
    if (numAccumulated < numAveraged)
    {
        if (!leftOut)
            log->info("Dynamic grid refinement: accumulated " + std::to_string(numAccumulated) + " of "
                      + std::to_string(numAveraged) + " iterations of the averaging window");
        return finish(false);
    }

    // ----- decision -----

    // evaluate all cells in parallel, recording for each cell and criterion whether the cell is a candidate or has
    // reached a limit, and for each cell the highest priority over the criteria for which it is a candidate;
    // each cell is evaluated by a single process, and the results are summed over the processes, so that they are
    // identical on all processes
    bool needsNeighbors = false;
    for (auto criterion : criteria) needsNeighbors |= criterion->needsNeighbors();
    Array priorityv(_numCells);
    Array flagv(static_cast<size_t>(_numCells) * numCriteria);
    find<ParallelFactory>()->parallelDistributed()->call(_numCells, [&](size_t firstIndex, size_t numIndices) {
        vector<int> neighbors;
        const vector<int> noNeighbors;
        for (size_t m = firstIndex; m != firstIndex + numIndices; ++m)
        {
            if (needsNeighbors) _tree->faceNeighbors(m, neighbors);
            int gainedLevels = _tree->cellLevel(m) - _tree->initialCellLevel(m);
            double size = _tree->cellSize(m);
            for (int c = 0; c != numCriteria; ++c)
            {
                auto criterion = criteria[c];

                // skip the cell if one of the conditions of the criterion does not hold
                bool holds = true;
                for (const auto& entry : conditionWindows[c])
                {
                    if (!entry.first->holds(entry.second->mean(m)))
                    {
                        holds = false;
                        break;
                    }
                }
                if (!holds) continue;

                // skip the cell if the criterion does not flag it
                double excess =
                    criterion->excess(m, *fieldWindows[c], criterion->needsNeighbors() ? neighbors : noNeighbors);
                if (!(excess > 0.)) continue;

                // check the limits of the criterion; if none has been reached, the cell is a candidate
                int flags = 0;
                int maxExtraLevels = criterion->maxExtraLevels();
                if (maxExtraLevels && gainedLevels >= maxExtraLevels) flags |= levelFlag;
                double minColumnDensity = criterion->minColumnDensity();
                double columnDensity = criterion->numberDensity(m) * size;
                if (minColumnDensity && !(columnDensity >= minColumnDensity)) flags |= columnDensityFlag;
                if (!flags)
                {
                    flags = candidateFlag;
                    double priority = excess * (minColumnDensity ? columnDensity / minColumnDensity : 1.);
                    priorityv[m] = max(priorityv[m], priority);
                }
                flagv[m * numCriteria + c] = flags;
            }
        }
    });
    ProcessManager::sumToAll(priorityv);
    ProcessManager::sumToAll(flagv);

    // count the candidates and the cells blocked by each limit, and collect the candidates in order of cell index
    vector<int> candidatev;
    vector<int> numCandidatesv(numCriteria), numLevelv(numCriteria), numColumnDensityv(numCriteria);
    for (int m = 0; m != _numCells; ++m)
    {
        bool isCandidate = false;
        for (int c = 0; c != numCriteria; ++c)
        {
            int flags = static_cast<int>(flagv[static_cast<size_t>(m) * numCriteria + c]);
            if (flags & candidateFlag)
            {
                numCandidatesv[c]++;
                isCandidate = true;
            }
            if (flags & levelFlag) numLevelv[c]++;
            if (flags & columnDensityFlag) numColumnDensityv[c]++;
        }
        if (isCandidate) candidatev.push_back(m);
    }
    int numCandidates = candidatev.size();
    candidates.set(numCandidates);
    log->info("Dynamic grid refinement: " + std::to_string(numCandidates) + " candidate cells out of "
              + std::to_string(_numCells));
    for (int c = 0; c != numCriteria; ++c)
    {
        criterionCandidates[c]->set(numCandidatesv[c]);
        blockedByLevel[c]->set(numLevelv[c]);
        blockedByColumnDensity[c]->set(numColumnDensityv[c]);
        log->info("  " + criterionDescription(criteria[c]) + ": " + std::to_string(numCandidatesv[c])
                  + " candidates; blocked by level limit: " + std::to_string(numLevelv[c])
                  + ", by column density limit: " + std::to_string(numColumnDensityv[c]));
    }

    // decide whether to start a refinement round, and how many candidates to subdivide
    int numSelected = 0;
    string reason;
    if (numCandidates)
    {
        int maxRounds = options->maxRounds();
        int maxCells = options->maxCells();
        if (maxRounds && numRounds >= maxRounds)
        {
            reason = "the maximum number of refinement rounds has been reached";
        }
        else if (maxIterations - iteration < numLeftOut + numAveraged)
        {
            reason = "too few iterations remain in the loop";
        }
        else
        {
            numSelected = numCandidates;
            if (maxCells)
            {
                int room = maxCells > _numCells ? (maxCells - _numCells) / (_tree->numChildren() - 1) : 0;
                if (room < numSelected) numSelected = room;
                if (!numSelected) reason = "the maximum number of cells has been reached";
            }
        }
    }

    // if the cell cap does not leave room for all candidates, select those with the highest priority; among
    // candidates with equal priority, the larger cells go first, and among those the cells with the lower index
    if (numSelected && numSelected < numCandidates)
    {
        std::sort(candidatev.begin(), candidatev.end(), [this, &priorityv](int m1, int m2) {
            if (priorityv[m1] != priorityv[m2]) return priorityv[m1] > priorityv[m2];
            double size1 = _tree->cellSize(m1);
            double size2 = _tree->cellSize(m2);
            if (size1 != size2) return size1 > size2;
            return m1 < m2;
        });
        lowest.set(priorityv[candidatev[numSelected - 1]]);
        candidatev.resize(numSelected);
        std::sort(candidatev.begin(), candidatev.end());
        reason = "the maximum number of cells has been reached";
    }

    // perform the refinement round, if any
    if (numSelected)
    {
        int numOldCells = _numCells;
        auto parentCells = _tree->subdivide(candidatev);
        growCellData(candidatev, parentCells);
        rounds.set(numRounds + 1);
        subdivided.set(numSelected);

        // record the aggregates again, so that they describe the refined grid
        recordAggregates();

        log->info("  Subdivided " + std::to_string(numSelected) + " cells, adding "
                  + std::to_string(_numCells - numOldCells) + " cells, for a total of " + std::to_string(_numCells)
                  + " cells (" + StringUtils::toString(numPackets / numOldCells, 'g', 3)
                  + " photon packets per cell in this iteration)");
    }
    if (numCandidates > numSelected)
        log->warning("Dynamic grid refinement left " + std::to_string(numCandidates - numSelected)
                     + " candidate cells unsplit because " + reason);

    // reset the windows, so that a new window starts with the next iteration that is not left out
    for (auto window : fieldWindows) window->reset();
    for (auto& entries : conditionWindows)
        for (auto& entry : entries) entry.second->reset();

    // the refinement has settled if the decision did not lead to a refinement round
    return finish(numSelected == 0);
}

////////////////////////////////////////////////////////////////////

double MediumSystem::callWithMaterialState(std::function<double(const MaterialState*)> callback, int m, int h) const
{
    MaterialState mst(_state, m, h);
    return callback(&mst);
}

////////////////////////////////////////////////////////////////////
