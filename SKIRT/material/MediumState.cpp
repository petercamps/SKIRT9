/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MediumState.hpp"
#include "FatalError.hpp"
#include "ProcessManager.hpp"
#include "StateVariableList.hpp"

//////////////////////////////////////////////////////////////////////

void MediumState::initConfiguration(int numCells, int numMedia)
{
    _numCells = numCells;
    _numMedia = numMedia;

    _off_dens.resize(_numMedia);
    _off_cust.resize(_numMedia);
}

//////////////////////////////////////////////////////////////////////

void MediumState::initCommonStateVariables(const vector<StateVariable>& variables)
{
    for (const StateVariable& variable : variables)
    {
        switch (variable.identifier())
        {
            case StateVariable::Identifier::Volume: _off_volu = _nextOffset++; break;
            case StateVariable::Identifier::BulkVelocity:
                _off_velo = _nextOffset;
                _nextOffset += 3;
                break;
            case StateVariable::Identifier::MagneticField:
                _off_mfld = _nextOffset;
                _nextOffset += 3;
                break;
            case StateVariable::Identifier::NumberDensity:
            case StateVariable::Identifier::Custom:
                throw FATALERROR("Requesting common state variable of unsupported type");
        }
    }
}

//////////////////////////////////////////////////////////////////////

void MediumState::initSpecificStateVariables(const StateVariableList& variables)
{
    // the list starts with the number density and assigns consecutive custom indices starting at zero,
    // so that the custom values of the component are stored contiguously after the number density
    for (const StateVariable& variable : variables)
    {
        switch (variable.identifier())
        {
            case StateVariable::Identifier::NumberDensity:
                _off_dens[_nextComponent] = _nextOffset++;
                _off_cust[_nextComponent] = _nextOffset;
                break;
            case StateVariable::Identifier::Custom:
                if (variable.kind() == QuantityKind::Extensive)
                    for (int k = 0; k != variable.count(); ++k) _off_extensive.push_back(_nextOffset + k);
                _nextOffset += variable.count();
                break;
            case StateVariable::Identifier::Volume:
            case StateVariable::Identifier::BulkVelocity:
            case StateVariable::Identifier::MagneticField:
                throw FATALERROR("Requesting specific state variable of unsupported type");
        }
    }
    _nextComponent++;
}

//////////////////////////////////////////////////////////////////////

size_t MediumState::initAllocate()
{
    if (_nextComponent != _numMedia) throw FATALERROR("Failed to request state variables for all medium components");
    _numVars = _nextOffset;

    size_t numAlloc = static_cast<size_t>(_numVars) * static_cast<size_t>(_numCells);
    _data.resize(numAlloc);
    return numAlloc;
}

//////////////////////////////////////////////////////////////////////

void MediumState::initCommunicate()
{
    ProcessManager::sumToAll(_data);
}

//////////////////////////////////////////////////////////////////////

std::pair<int, int> MediumState::synchronize(const vector<UpdateStatus>& cellFlags)
{
    int numUpdated = 0;
    int numNotConverged = 0;

    if (ProcessManager::isMultiProc())
    {
        auto producer = [this, &cellFlags, &numUpdated, &numNotConverged](vector<double>& data) {
            for (int m = 0; m != _numCells; ++m)
            {
                if (cellFlags[m].isUpdated())
                {
                    // cell index
                    data.push_back(m);
                    // state variables
                    double* first = &_data[_numVars * m];
                    data.insert(data.end(), first, first + _numVars);
                    // update status
                    numUpdated++;
                    if (cellFlags[m].isConverged())
                    {
                        data.push_back(0.);
                    }
                    else
                    {
                        data.push_back(1.);
                        numNotConverged++;
                    }
                }
            }
        };
        auto consumer = [this, &numUpdated, &numNotConverged](const vector<double>& data) {
            for (auto in = data.begin(); in != data.end(); in += (2 + _numVars))
            {
                // cell index
                int m = *in;
                // state variables
                std::copy(in + 1, in + 1 + _numVars, &_data[_numVars * m]);
                // update status
                numUpdated++;
                if (*(in + 1 + _numVars)) numNotConverged++;
            }
        };
        ProcessManager::broadcastAllToAll(producer, consumer);
    }
    else
    {
        for (int m = 0; m != _numCells; ++m)
        {
            if (cellFlags[m].isUpdated()) numUpdated++;
            if (!cellFlags[m].isConverged()) numNotConverged++;
        }
    }
    return std::make_pair(numUpdated, numNotConverged);
}

//////////////////////////////////////////////////////////////////////

void MediumState::appendCells(const vector<int>& parentCells)
{
    // copy all variables of the parent cell to each new cell; the parents precede the new cells, so that copying
    // after growing the vector never reads values of a new cell
    size_t numOld = _numCells;
    size_t numNew = parentCells.size();
    _data.resize(_numVars * (numOld + numNew));
    for (size_t i = 0; i != numNew; ++i)
    {
        auto parent = _data.begin() + _numVars * parentCells[i];
        std::copy(parent, parent + _numVars, _data.begin() + _numVars * (numOld + i));
    }
    _numCells = numOld + numNew;
}

//////////////////////////////////////////////////////////////////////

void MediumState::resizeCell(int m, double volume)
{
    double oldVolume = this->volume(m);
    double fraction = oldVolume > 0. ? volume / oldVolume : 0.;
    for (int offset : _off_extensive) _data[_numVars * m + offset] *= fraction;
    setVolume(m, volume);
}

//////////////////////////////////////////////////////////////////////

void MediumState::reserveCells(int numCells)
{
    _data.reserve(_numVars * static_cast<size_t>(numCells));
}

//////////////////////////////////////////////////////////////////////

vector<double> MediumState::volumeIntegrals(const vector<std::pair<int, int>>& variables) const
{
    vector<double> result(variables.size(), 0.);
    for (int m = 0; m != _numCells; ++m)
    {
        double volume = _data[_numVars * m + _off_volu];
        for (size_t k = 0; k != variables.size(); ++k)
            result[k] += _data[_numVars * m + _off_cust[variables[k].first] + variables[k].second] * volume;
    }
    return result;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setVolume(int m, double value)
{
    _data[_numVars * m + _off_volu] = value;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setBulkVelocity(int m, Vec value)
{
    size_t i = _numVars * m + _off_velo;
    _data[i] = value.x();
    _data[i + 1] = value.y();
    _data[i + 2] = value.z();
}

//////////////////////////////////////////////////////////////////////

void MediumState::setMagneticField(int m, Vec value)
{
    size_t i = _numVars * m + _off_mfld;
    _data[i] = value.x();
    _data[i + 1] = value.y();
    _data[i + 2] = value.z();
}

//////////////////////////////////////////////////////////////////////

void MediumState::setNumberDensity(int m, int h, double value)
{
    _data[_numVars * m + _off_dens[h]] = value;
}

//////////////////////////////////////////////////////////////////////

void MediumState::setCustom(int m, int h, int i, double value)
{
    _data[_numVars * m + _off_cust[h] + i] = value;
}

//////////////////////////////////////////////////////////////////////
