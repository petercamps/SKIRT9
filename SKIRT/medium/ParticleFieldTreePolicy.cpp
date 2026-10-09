/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ParticleFieldTreePolicy.hpp"
#include "FatalError.hpp"
#include "ParticleSnapshot.hpp"
#include "TreeNodeEvaluation.hpp"

////////////////////////////////////////////////////////////////////

ParticleFieldTreePolicy::~ParticleFieldTreePolicy()
{
    delete _snapshot;
}

////////////////////////////////////////////////////////////////////

void ParticleFieldTreePolicy::setupSelfAfter()
{
    TreePolicy::setupSelfAfter();

    // import the particles, with a mass density policy that keeps the masses as they are,
    // so that the snapshot offers the total mass, the mass in a box, and the density at a position
    _snapshot = new ParticleSnapshot;
    _snapshot->open(this, filename(), "smoothed particles");
    _snapshot->importPosition();
    _snapshot->importSize();
    _snapshot->importMass();
    _snapshot->setSmoothingKernel(smoothingKernel());
    _snapshot->setMassDensityPolicy(1., 0., false);
    _snapshot->readAndClose();

    if (!(_snapshot->mass() > 0.)) throw FATALERROR("The total mass of the smoothed particle field is not positive");
    _hasMassInBox = smoothingKernel()->hasMassInBox();
    _maxMass = maxFraction() * _snapshot->mass();
}

////////////////////////////////////////////////////////////////////

bool ParticleFieldTreePolicy::needsSubdivide(TreeNodeEvaluation& node) const
{
    double mass = 0.;
    if (_hasMassInBox)
    {
        mass = _snapshot->massInBox(node.box());
    }
    else
    {
        const auto& positionv = node.samplePositions();
        for (const auto& position : positionv) mass += _snapshot->density(position);
        mass *= node.box().volume() / positionv.size();
    }
    return mass > _maxMass;
}

////////////////////////////////////////////////////////////////////
