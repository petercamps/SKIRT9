/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParticleSnapshot.hpp"
#include "Log.hpp"
#include "NR.hpp"
#include "Random.hpp"
#include "SmoothedParticleGrid.hpp"
#include "SmoothingKernel.hpp"
#include "StringUtils.hpp"
#include "TextInFile.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

ParticleSnapshot::~ParticleSnapshot()
{
    delete _grid;
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::open(const SimulationItem* item, string filename, string description)
{
    Snapshot::open(item, filename, description);
    importPosition();
    importSize();
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::readAndClose()
{
    // read the particle info into memory
    // if the user configured a temperature cutoff, we need to skip the "hot" particles
    int numIgnored = 0;
    if (!hasTemperatureCutoff()) _propv = infile()->readAllRows();
    else
    {
        Array row;
        while (infile()->readRow(row))
        {
            if (row[temperatureIndex()] > maxTemperature()) numIgnored++;
            else _propv.push_back(row);
        }
    }

    // close the file
    Snapshot::readAndClose();

    // log the number of particles
    if (!numIgnored)
    {
        log()->info("  Number of particles: " + std::to_string(_propv.size()));
    }
    else
    {
        log()->info("  Number of particles ignored: " + std::to_string(numIgnored));
        log()->info("  Number of particles retained: " + std::to_string(_propv.size()));
    }

    // we can calculate mass and densities only if a policy has been set
    if (!hasMassDensityPolicy()) return;

    // build a list of compact smoothed particle objects that we can organize in a grid
    double totalOriginalMass = 0;
    double totalMetallicMass = 0;
    double totalEffectiveMass = 0;
    _pv.reserve(_propv.size());
    for (const Array& prop : _propv)
    {
        double originalMass = prop[massIndex()];
        double metallicMass = originalMass * (metallicityIndex()>=0 ? prop[metallicityIndex()] : 1.);
        double effectiveMass = metallicMass * multiplier();

        _pv.emplace_back(prop[positionIndex()+0], prop[positionIndex()+1], prop[positionIndex()+2],
                prop[sizeIndex()], effectiveMass);

        totalOriginalMass += originalMass;
        totalMetallicMass += metallicMass;
        totalEffectiveMass += effectiveMass;
    }

    // log mass statistics
    log()->info("  Total original mass : " +
                StringUtils::toString(units()->omass(totalOriginalMass)) + units()->umass());
    log()->info("  Total metallic mass : "
                + StringUtils::toString(units()->omass(totalMetallicMass)) + units()->umass());
    log()->info("  Total effective mass: "
                + StringUtils::toString(units()->omass(totalEffectiveMass)) + units()->umass());

    // if one of the total masses is negative, suppress the complete mass distribution
    if (totalOriginalMass < 0 || totalMetallicMass < 0 || totalEffectiveMass < 0)
    {
        log()->warning("  Total imported mass is negative; suppressing the complete mass distribition");
        _propv.clear();
        _pv.clear();
        return;         // abort
    }

    // remember the effective mass
    _mass = totalEffectiveMass;

    // if there are no particles, do not build the special structures for optimizing operations
    if (_pv.empty()) return;

    // construct a 3D-grid over the particle space, and create a list of particles that overlap each grid cell
    int gridsize = max(20, static_cast<int>(  pow(_pv.size(),1./3.)/5 ));
    string size = std::to_string(gridsize);
    log()->info("Constructing intermediate " + size + "x" + size + "x" + size + " grid for particles...");
    _grid = new SmoothedParticleGrid(_pv, gridsize);
    log()->info("  Smallest number of particles per cell: " + std::to_string(_grid->minParticlesPerCell()));
    log()->info("  Largest  number of particles per cell: " + std::to_string(_grid->maxParticlesPerCell()));
    log()->info("  Average  number of particles per cell: "
                + StringUtils::toString(_grid->totalParticles() / double(gridsize*gridsize*gridsize),'f',1));

    // construct a vector with the normalized cumulative particle densities
    NR::cdf(_cumrhov, _pv.size(), [this](int i){return _pv[i].mass();} );
}

////////////////////////////////////////////////////////////////////

void ParticleSnapshot::setSmoothingKernel(const SmoothingKernel* kernel)
{
    _kernel = kernel;
}

////////////////////////////////////////////////////////////////////

Box ParticleSnapshot::extent() const
{
    // if there are no particles, return an empty box
    if (_propv.empty()) return Box();

    // if there is a particle grid, ask it to return the extent (it is already calculated)
    if (_grid) return _grid->extent();

    // otherwise find the spatial range of the particles assuming a finite support kernel
    double xmin = + std::numeric_limits<double>::infinity();
    double xmax = - std::numeric_limits<double>::infinity();
    double ymin = + std::numeric_limits<double>::infinity();
    double ymax = - std::numeric_limits<double>::infinity();
    double zmin = + std::numeric_limits<double>::infinity();
    double zmax = - std::numeric_limits<double>::infinity();
    for (const Array& prop : _propv)
    {
        xmin = min(xmin, prop[positionIndex()+0] - prop[sizeIndex()]);
        xmax = max(xmax, prop[positionIndex()+0] + prop[sizeIndex()]);
        ymin = min(ymin, prop[positionIndex()+1] - prop[sizeIndex()]);
        ymax = max(ymax, prop[positionIndex()+1] + prop[sizeIndex()]);
        zmin = min(zmin, prop[positionIndex()+2] - prop[sizeIndex()]);
        zmax = max(zmax, prop[positionIndex()+2] + prop[sizeIndex()]);
    }
    return Box(xmin, ymin, zmin, xmax, ymax, zmax);
}

////////////////////////////////////////////////////////////////////

int ParticleSnapshot::numEntities() const
{
    return _propv.size();
}

////////////////////////////////////////////////////////////////////

Position ParticleSnapshot::position(int m) const
{
    return Position(_propv[m][positionIndex()+0], _propv[m][positionIndex()+1], _propv[m][positionIndex()+2]);
}

////////////////////////////////////////////////////////////////////

Vec ParticleSnapshot::velocity(int m) const
{
    return Vec(_propv[m][velocityIndex()+0], _propv[m][velocityIndex()+1], _propv[m][velocityIndex()+2]);
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::density(Position bfr) const
{
    double sum = 0.;
    if (_grid) for (const SmoothedParticle* p : _grid->particlesFor(bfr))
    {
        double u = (bfr - p->center()).norm() / p->radius();
        sum += _kernel->density(u) * p->mass();
    }
    return sum > 0. ? sum : 0.;     // guard against negative densities
}

////////////////////////////////////////////////////////////////////

double ParticleSnapshot::mass() const
{
    return _mass;
}

////////////////////////////////////////////////////////////////////

Position ParticleSnapshot::generatePosition(int m) const
{
    // get center position and size for this particle
    Position rc(_propv[m][positionIndex()+0], _propv[m][positionIndex()+1], _propv[m][positionIndex()+2]);
    double h = _propv[m][sizeIndex()];

    // sample random position inside the smoothed unit volume
    double u = _kernel->generateRadius();
    Direction k = random()->direction();

    return Position(rc + k*u*h);
}

////////////////////////////////////////////////////////////////////

Position ParticleSnapshot::generatePosition() const
{
    // if there are no particles, return the origin
    if (_propv.empty()) return Position();

    // select a particle according to its mass contribution
    int m = NR::locateClip(_cumrhov, random()->uniform());

    return generatePosition(m);
}

////////////////////////////////////////////////////////////////////