/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef QUANTITYKIND_HPP
#define QUANTITYKIND_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** The QuantityKind enumeration lists the kinds of a quantity that is imported from a snapshot,
    as described by a SnapshotParameter object. The kind determines how the values of several
    entities (e.g., overlapping smoothed particles) at a given position, and of several samples in
    a spatial cell, are combined:

    - \c Intensive: a continuous quantity that does not depend on the amount of material, such as
      a temperature or an abundance ratio. The values of overlapping entities are averaged,
      weighted by the kernel weight times the mass of each entity, and the values of the samples
      in a cell are averaged, weighted by density.

    - \c Categorical: a label or an index, for which an average has no meaning, such as the index
      of a dust mix. The value is taken from the overlapping entity with the largest kernel
      weight, and from the sample in a cell with the highest density.

    - \c Extensive: a quantity that is proportional to the amount of material, such as a mass.
      Combining extensive quantities is not supported.

    A quantity imported for a material mix or a material mix family is intensive or categorical; a
    material mix that needs an extensive quantity requests the corresponding intensive quantity
    instead, such as a number density rather than a number of particles. */
enum class QuantityKind { Intensive, Categorical, Extensive };

////////////////////////////////////////////////////////////////////

#endif
