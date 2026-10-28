/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SECONDARYSOURCE_HPP
#define SECONDARYSOURCE_HPP

#include "SimulationItem.hpp"
class PhotonPacket;

//////////////////////////////////////////////////////////////////////

/** SecondarySource is the abstract base class for helper classes that handle launching photon
    packets from a secondary source of a given type, i.e. dust or gas. The class inherits
    SimulationItem so that, for example, it can easily use the find() template function, but it is
    not part of the configurable simulation item hierarchy. Instead, the SecondarySourceSystem
    class creates and holds SecondarySource subclass instances as required depending on the
    simulation's configuration.

    For more information on the operation of this class and its subclasses, see the
    SecondarySourceSystem class. */
class SecondarySource : public SimulationItem
{
    //============= Construction - Setup - Destruction =============

public:
    /** This constructor creates a SecondarySource subclass instance. Before the constructor
        returns, the newly created object is hooked up as a child to the specified parent in the
        simulation hierarchy (so it will automatically be deleted). */
    explicit SecondarySource(SimulationItem* parent);

    //======================== Other Functions =======================

public:
    /** This function calculates and stores the bolometric luminosities for this source in each
        spatial cell of the simulation, and returns the grand-total bolometric luminosity for this
        source. */
    virtual double prepareLuminosities() = 0;

    /** This function prepares the mapping of history indices to sources, given the range of
        history indices allocated to this source . */
    virtual void preparePacketMap(size_t firstIndex, size_t numIndices) = 0;

    /** This function causes the photon packet \em pp to be launched for this source from one of
        the cells in the spatial grid using the given history index and luminosity. */
    virtual void launch(PhotonPacket* pp, size_t historyIndex, double L) const = 0;

    //======================== Launch preparations =======================

protected:
    /** This function starts a new launch preparation for this source, and must be called at the
        start of the prepareLuminosities() function of each subclass. It assigns the source a
        preparation number that differs from the number of any earlier preparation of any secondary
        source in the process. The per-thread caches that a subclass uses while launching photon
        packets compare the preparation number with the one for which they hold information, so
        that they never use information calculated for another source, or for an earlier
        preparation, for example before the number of cells changed through dynamic grid refinement
        or before the radiation field was updated in a new iteration. */
    void beginPreparation();

    /** This function returns the preparation number assigned by the most recent call to
        beginPreparation(). */
    size_t preparation() const { return _preparation; }

    //======================== Data Members ========================

private:
    size_t _preparation{0};
};

////////////////////////////////////////////////////////////////

#endif
