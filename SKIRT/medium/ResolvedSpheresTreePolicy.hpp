/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef RESOLVEDSPHERESTREEPOLICY_HPP
#define RESOLVEDSPHERESTREEPOLICY_HPP

#include "BoxSearch.hpp"
#include "TreePolicy.hpp"

//////////////////////////////////////////////////////////////////////

/** ResolvedSpheresTreePolicy is a tree subdivision policy that ensures that the grid resolves
    each sphere in a list of spheres loaded from a text column file. Specifically, the policy
    subdivides a node if it intersects one of the spheres, enlarged by a factor \em reach, and if
    the node is wider than the radius of that sphere divided by \em numBins in any of the spatial
    directions. As a result, the grid has at least \em numBins cells across the radius of each
    sphere in each spatial direction, out to a distance of \em reach radii from the sphere's
    center, unless the maximum level of the grid is reached first.

    The input file contains a sphere on each line, with columns for the \f$x\f$, \f$y\f$, and
    \f$z\f$ coordinates of its center and for its radius. If the \em importNumBins flag is
    enabled, the next column specifies the number of bins for the sphere, and if the \em
    importReach flag is enabled, the next column specifies its reach. A positive imported value
    overrides the value configured for the policy for that sphere, while a zero or negative value
    selects the configured value. This allows specifying a different resolution or reach for some
    of the spheres. The default units are pc for the coordinates and the radius; this can be
    overridden by providing column information in the header of the file, as described for the
    TextInFile class. The number of bins and the reach are dimensionless. Spheres with a radius
    that is not positive are ignored.

    During setup, the policy organizes the enlarged spheres in a search structure (see the
    BoxSearch class), so that the spheres intersecting a given node can be found efficiently. */
class ResolvedSpheresTreePolicy : public TreePolicy
{
    ITEM_CONCRETE(ResolvedSpheresTreePolicy, TreePolicy,
                  "a tree subdivision policy resolving a list of spheres loaded from file")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ResolvedSpheresTreePolicy, "Level2")

        PROPERTY_STRING(filename, "the name of the file with the spheres")

        PROPERTY_INT(numBins, "the minimum number of cells across the radius of each sphere")
        ATTRIBUTE_MIN_VALUE(numBins, "1")
        ATTRIBUTE_MAX_VALUE(numBins, "10000")
        ATTRIBUTE_DEFAULT_VALUE(numBins, "10")

        PROPERTY_DOUBLE(reach, "the distance from the center, in sphere radii, out to which to resolve each sphere")
        ATTRIBUTE_MIN_VALUE(reach, "]0")
        ATTRIBUTE_MAX_VALUE(reach, "100]")
        ATTRIBUTE_DEFAULT_VALUE(reach, "1")

        PROPERTY_BOOL(importNumBins, "import a column with the number of bins for each sphere")
        ATTRIBUTE_DEFAULT_VALUE(importNumBins, "false")
        ATTRIBUTE_DISPLAYED_IF(importNumBins, "Level3")

        PROPERTY_BOOL(importReach, "import a column with the reach for each sphere")
        ATTRIBUTE_DEFAULT_VALUE(importReach, "false")
        ATTRIBUTE_DISPLAYED_IF(importReach, "Level3")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function imports the spheres and organizes them in a search structure. */
    void setupSelfBefore() override;

    //======================== Other Functions =======================

public:
    /** This function returns true if the specified node intersects one of the enlarged spheres
        and is wider than the required cell size for that sphere in any spatial direction, and false
        otherwise. */
    bool needsSubdivide(TreeNodeEvaluation& node) const override;

    //======================== Data Members ========================

private:
    // data members initialized by setupSelfBefore()
    vector<Vec> _centerv;    // the center of each sphere
    vector<double> _reachv;  // the radius of each sphere multiplied by its reach
    vector<double> _sizev;   // the maximum cell size for each sphere, i.e. its radius divided by its number of bins
    BoxSearch _search;       // the search structure for the enlarged spheres
};

//////////////////////////////////////////////////////////////////////

#endif
