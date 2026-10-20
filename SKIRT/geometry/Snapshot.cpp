/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "Snapshot.hpp"
#include "ColumnInFile.hpp"
#include "EntityCollection.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

Snapshot::Snapshot() {}

////////////////////////////////////////////////////////////////////

Snapshot::~Snapshot() {}

////////////////////////////////////////////////////////////////////

void Snapshot::open(const SimulationItem* item, const string& filename, const string& description)
{
    _infile = createInFile(item, filename, description);
    setContext(item);
}

////////////////////////////////////////////////////////////////////

std::unique_ptr<ColumnInFile> Snapshot::createInFile(const SimulationItem* item, const string& filename,
                                                     const string& description)
{
    return std::make_unique<ColumnInFile>(item, filename, description);
}

////////////////////////////////////////////////////////////////////

void Snapshot::close()
{
    _infile.reset();
}

////////////////////////////////////////////////////////////////////

void Snapshot::readAndClose()
{
    close();
}

////////////////////////////////////////////////////////////////////

void Snapshot::setContext(const SimulationItem* item)
{
    _log = item->find<Log>();
    _units = item->find<Units>();
    _random = item->find<Random>();
}

////////////////////////////////////////////////////////////////////

void Snapshot::useColumns(const string& columns)
{
    _infile->useColumns(columns);
}

////////////////////////////////////////////////////////////////////

void Snapshot::setCoordinateSystem(CoordinateSystem coordinateSystem)
{
    _coordinateSystem = coordinateSystem;
}

////////////////////////////////////////////////////////////////////

void Snapshot::importPosition()
{
    _positionIndex = _nextIndex;
    _nextIndex += 3;
    switch (_coordinateSystem)
    {
        case CoordinateSystem::CARTESIAN:
            _infile->addColumn("position x", "length", "pc");
            _infile->addColumn("position y", "length", "pc");
            _infile->addColumn("position z", "length", "pc");
            break;
        case CoordinateSystem::CYLINDRICAL:
            _infile->addColumn("position R", "length", "pc");
            _infile->addColumn("position phi", "posangle", "deg");
            _infile->addColumn("position z", "length", "pc");
            break;
        case CoordinateSystem::SPHERICAL:
            _infile->addColumn("position r", "length", "pc");
            _infile->addColumn("position theta", "posangle", "deg");
            _infile->addColumn("position phi", "posangle", "deg");
            break;
    }
}

////////////////////////////////////////////////////////////////////

void Snapshot::importSize()
{
    _sizeIndex = _nextIndex++;
    _infile->addColumn("size h", "length", "pc");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importBox()
{
    _boxIndex = _nextIndex;
    _nextIndex += 6;
    switch (_coordinateSystem)
    {
        case CoordinateSystem::CARTESIAN:
            _infile->addColumn("box xmin", "length", "pc");
            _infile->addColumn("box ymin", "length", "pc");
            _infile->addColumn("box zmin", "length", "pc");
            _infile->addColumn("box xmax", "length", "pc");
            _infile->addColumn("box ymax", "length", "pc");
            _infile->addColumn("box zmax", "length", "pc");
            break;
        case CoordinateSystem::CYLINDRICAL:
            _infile->addColumn("box Rmin", "length", "pc");
            _infile->addColumn("box phimin", "posangle", "deg");
            _infile->addColumn("box zmin", "length", "pc");
            _infile->addColumn("box Rmax", "length", "pc");
            _infile->addColumn("box phimax", "posangle", "deg");
            _infile->addColumn("box zmax", "length", "pc");
            break;
        case CoordinateSystem::SPHERICAL:
            _infile->addColumn("box rmin", "length", "pc");
            _infile->addColumn("box thetamin", "posangle", "deg");
            _infile->addColumn("box phimin", "posangle", "deg");
            _infile->addColumn("box rmax", "length", "pc");
            _infile->addColumn("box thetamax", "posangle", "deg");
            _infile->addColumn("box phimax", "posangle", "deg");
            break;
    }
}

////////////////////////////////////////////////////////////////////

void Snapshot::importMassDensity()
{
    _densityIndex = _nextIndex++;
    _infile->addColumn("mass density", "massvolumedensity", "Msun/pc3");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importMass()
{
    _massIndex = _nextIndex++;
    _infile->addColumn("mass", "mass", "Msun");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importNumberDensity()
{
    _densityIndex = _nextIndex++;
    _infile->addColumn("number density", "numbervolumedensity", "1/cm3");
    _holdsNumber = true;
}

////////////////////////////////////////////////////////////////////

void Snapshot::importNumber()
{
    _massIndex = _nextIndex++;
    _infile->addColumn("number");
    _holdsNumber = true;
}

////////////////////////////////////////////////////////////////////

void Snapshot::importCurrentMass()
{
    _currentMassIndex = _nextIndex++;
    _infile->addColumn("current mass", "mass", "Msun");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importVelocity()
{
    _velocityIndex = _nextIndex;
    _nextIndex += 3;
    switch (_coordinateSystem)
    {
        case CoordinateSystem::CARTESIAN:
            _infile->addColumn("velocity x", "velocity", "km/s");
            _infile->addColumn("velocity y", "velocity", "km/s");
            _infile->addColumn("velocity z", "velocity", "km/s");
            break;
        case CoordinateSystem::CYLINDRICAL:
            _infile->addColumn("velocity R", "velocity", "km/s");
            _infile->addColumn("velocity phi", "velocity", "km/s");
            _infile->addColumn("velocity z", "velocity", "km/s");
            break;
        case CoordinateSystem::SPHERICAL:
            _infile->addColumn("velocity r", "velocity", "km/s");
            _infile->addColumn("velocity theta", "velocity", "km/s");
            _infile->addColumn("velocity phi", "velocity", "km/s");
            break;
    }
}

////////////////////////////////////////////////////////////////////

void Snapshot::importVelocityDispersion()
{
    _velocityDispersionIndex = _nextIndex++;
    _infile->addColumn("velocity dispersion", "velocity", "km/s");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importMagneticField()
{
    _magneticFieldIndex = _nextIndex;
    _nextIndex += 3;
    switch (_coordinateSystem)
    {
        case CoordinateSystem::CARTESIAN:
            _infile->addColumn("magnetic field x", "magneticfield", "uG");
            _infile->addColumn("magnetic field y", "magneticfield", "uG");
            _infile->addColumn("magnetic field z", "magneticfield", "uG");
            break;
        case CoordinateSystem::CYLINDRICAL:
            _infile->addColumn("magnetic field R", "magneticfield", "uG");
            _infile->addColumn("magnetic field phi", "magneticfield", "uG");
            _infile->addColumn("magnetic field z", "magneticfield", "uG");
            break;
        case CoordinateSystem::SPHERICAL:
            _infile->addColumn("magnetic field r", "magneticfield", "uG");
            _infile->addColumn("magnetic field theta", "magneticfield", "uG");
            _infile->addColumn("magnetic field phi", "magneticfield", "uG");
            break;
    }
}

////////////////////////////////////////////////////////////////////

void Snapshot::importBias()
{
    _biasIndex = _nextIndex++;
    _infile->addColumn("bias");
}

////////////////////////////////////////////////////////////////////

void Snapshot::importParameters(const vector<SnapshotParameter>& parameters)
{
    _parametersIndex = _nextIndex;
    _numParameters = parameters.size();
    for (const auto& param : parameters)
    {
        _infile->addColumn(param.description(), param.quantity(), param.defaultUnit());
        _parameterKinds.push_back(param.kind());

        // discover standard parameter types and set the corresponding index so these values can be probed
        switch (param.identifier())
        {
            case SnapshotParameter::Identifier::InitialMass: _initialMassIndex = _nextIndex; break;
            case SnapshotParameter::Identifier::CurrentMass: _currentMassIndex = _nextIndex; break;
            case SnapshotParameter::Identifier::Metallicity: _metallicityIndex = _nextIndex; break;
            case SnapshotParameter::Identifier::Age: _ageIndex = _nextIndex; break;
            case SnapshotParameter::Identifier::Temperature: _temperatureIndex = _nextIndex; break;
            case SnapshotParameter::Identifier::Custom: break;
        }
        _nextIndex++;
    }
}

////////////////////////////////////////////////////////////////////

void Snapshot::setMassDensityPolicy(double multiplier)
{
    _multiplier = multiplier;
    _hasDensityPolicy = true;
}

////////////////////////////////////////////////////////////////////

void Snapshot::setNeedGetEntities()
{
    _needGetEntities = true;
}

////////////////////////////////////////////////////////////////////

void Snapshot::logMassStatistics(double totalOriginalMass, double totalEffectiveMass)
{
    auto massUnit = holdsNumber() ? "" : units()->umass();
    auto massToString = [this](double mass) {
        return StringUtils::toString(holdsNumber() ? mass : units()->omass(mass), 'e', 4);
    };
    log()->info("  Total original mass : " + massToString(totalOriginalMass) + " " + massUnit);
    log()->info("  Total effective mass: " + massToString(totalEffectiveMass) + " " + massUnit);
}

////////////////////////////////////////////////////////////////////

void Snapshot::calculateDensityAndMass(Array& rhov, Array& cumrhov, double& mass)
{
    // resize the density array and allocate a temporary mass vector
    int numCells = numEntities();
    rhov.resize(numCells);
    Array Mv(numCells);

    // initialize statistics
    double totalOriginalMass = 0.;
    double totalEffectiveMass = 0.;

    // loop over all cells
    for (int m = 0; m != numCells; ++m)
    {
        const Array& prop = properties(m);

        // original mass is zero if imported mass/density is not positive;
        // use density or mass or both, depending on availability
        double V = volume(m);
        double originalDensity = max(0., densityIndex() >= 0 ? prop[densityIndex()] : prop[massIndex()] / V);
        double originalMass = max(0., massIndex() >= 0 ? prop[massIndex()] : prop[densityIndex()] * V);

        // determine effective mass
        double effectiveDensity = originalDensity * multiplier();
        double effectiveMass = originalMass * multiplier();

        // store density and mass for this cell
        rhov[m] = effectiveDensity;
        Mv[m] = effectiveMass;

        // accumulate statistics
        totalOriginalMass += originalMass;
        totalEffectiveMass += effectiveMass;
    }

    // construct and store a vector with the normalized cumulative site densities
    if (numCells) NR::cdf(cumrhov, Mv);

    // store the total effective mass
    mass = totalEffectiveMass;

    // log mass statistics
    logMassStatistics(totalOriginalMass, totalEffectiveMass);
}

////////////////////////////////////////////////////////////////////

double Snapshot::massForAveraging(int /*m*/) const
{
    return 1.;
}

////////////////////////////////////////////////////////////////////

double Snapshot::volume() const
{
    return extent().volume();
}

////////////////////////////////////////////////////////////////////

double Snapshot::initialMass(int m) const
{
    return properties(m)[initialMassIndex()];
}

////////////////////////////////////////////////////////////////////

double Snapshot::currentMass(int m) const
{
    return properties(m)[currentMassIndex()];
}

////////////////////////////////////////////////////////////////////

double Snapshot::metallicity(int m) const
{
    return properties(m)[metallicityIndex()];
}

////////////////////////////////////////////////////////////////////

double Snapshot::age(int m) const
{
    return properties(m)[ageIndex()];
}

////////////////////////////////////////////////////////////////////

double Snapshot::temperature(int m) const
{
    return properties(m)[temperatureIndex()];
}

////////////////////////////////////////////////////////////////////

Vec Snapshot::velocity(int m) const
{
    const auto& propv = properties(m);
    return Vec(propv[velocityIndex() + 0], propv[velocityIndex() + 1], propv[velocityIndex() + 2]);
}

////////////////////////////////////////////////////////////////////

Vec Snapshot::velocity(Position bfr) const
{
    thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
    getEntities(entities, bfr);
    return entities.averageValue([this](int m) { return velocity(m); }, [this](int m) { return massForAveraging(m); });
}

////////////////////////////////////////////////////////////////////

double Snapshot::velocityDispersion(int m) const
{
    return properties(m)[velocityDispersionIndex()];
}

////////////////////////////////////////////////////////////////////

Vec Snapshot::magneticField(int m) const
{
    const auto& propv = properties(m);
    return Vec(propv[magneticFieldIndex() + 0], propv[magneticFieldIndex() + 1], propv[magneticFieldIndex() + 2]);
}

////////////////////////////////////////////////////////////////////

Vec Snapshot::magneticField(Position bfr) const
{
    thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
    getEntities(entities, bfr);
    return entities.averageValue([this](int m) { return magneticField(m); },
                                 [this](int m) { return massForAveraging(m); });
}

////////////////////////////////////////////////////////////////////

double Snapshot::bias(int m) const
{
    return properties(m)[biasIndex()];
}

////////////////////////////////////////////////////////////////////

void Snapshot::parameters(int m, Array& params) const
{
    int n = numParameters();
    params.resize(n);
    const auto& propv = properties(m);
    for (int i = 0; i != n; ++i) params[i] = propv[parametersIndex() + i];
}

////////////////////////////////////////////////////////////////////

void Snapshot::parameters(Position bfr, Array& params) const
{
    thread_local EntityCollection entities;  // can be reused for all queries in a given execution thread
    getEntities(entities, bfr);

    // without entities, return zero values; for a single entity, return its values
    int n = numParameters();
    auto first = entities.begin();
    if (first == entities.end())
    {
        params.resize(n);
        return;
    }
    if (std::next(first) == entities.end())
    {
        parameters(first->first, params);
        return;
    }

    // in a single pass over the entities, accumulate the weighted values of the intensive parameters,
    // and locate the entity with the largest weight for the categorical parameters
    params.resize(n);
    double sumw = 0.;
    double wmax = 0.;
    int mmax = -1;
    for (const auto& [m, wm] : entities)
    {
        double w = massForAveraging(m) * wm;
        const auto& propv = properties(m);
        for (int i = 0; i != n; ++i)
            if (_parameterKinds[i] == QuantityKind::Intensive) params[i] += propv[parametersIndex() + i] * w;
        sumw += w;
        if (wm > wmax)
        {
            wmax = wm;
            mmax = m;
        }
    }

    // finalize the averages, guarding against weights that add up to zero, and copy the categorical values;
    // an entity always has a positive weight in the collection, so that the entity with the largest weight exists
    const auto& propmax = properties(mmax);
    for (int i = 0; i != n; ++i)
    {
        if (_parameterKinds[i] == QuantityKind::Intensive)
            params[i] = sumw > 0. ? params[i] / sumw : 0.;
        else
            params[i] = propmax[parametersIndex() + i];
    }
}

////////////////////////////////////////////////////////////////////

namespace
{
    // the number of samples used for integrating the surface densities
    const int NSAMPLES = 10000;
}

////////////////////////////////////////////////////////////////////

double Snapshot::SigmaX() const
{
    // determine a small value relative to the domain extent;
    // we integrate along a small offset from the axis to avoid cell borders
    double eps = 1e-12 * extent().diagonal();

    double sum = 0;
    double xmin = extent().xmin();
    double xmax = extent().xmax();
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(xmin + k * (xmax - xmin) / NSAMPLES, eps, eps));
    }
    return (sum / NSAMPLES) * (xmax - xmin);
}

////////////////////////////////////////////////////////////////////

double Snapshot::SigmaY() const
{
    // determine a small value relative to the domain extent;
    // we integrate along a small offset from the axis to avoid cell borders
    double eps = 1e-12 * extent().diagonal();

    double sum = 0;
    double ymin = extent().ymin();
    double ymax = extent().ymax();
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(eps, ymin + k * (ymax - ymin) / NSAMPLES, eps));
    }
    return (sum / NSAMPLES) * (ymax - ymin);
}

////////////////////////////////////////////////////////////////////

double Snapshot::SigmaZ() const
{
    // determine a small value relative to the domain extent;
    // we integrate along a small offset from the axis to avoid cell borders
    double eps = 1e-12 * extent().diagonal();

    double sum = 0;
    double zmin = extent().zmin();
    double zmax = extent().zmax();
    for (int k = 0; k < NSAMPLES; k++)
    {
        sum += density(Position(eps, eps, zmin + k * (zmax - zmin) / NSAMPLES));
    }
    return (sum / NSAMPLES) * (zmax - zmin);
}

////////////////////////////////////////////////////////////////////
