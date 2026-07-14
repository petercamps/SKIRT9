/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef CLUMPYTORUSSPATIALGRID_HPP
#define CLUMPYTORUSSPATIALGRID_HPP

#include "SpatialGrid.hpp"

////////////////////////////////////////////////////////////////////

/** An instance of the ClumpyTorusSpatialGrid class represents a specialty grid that handles a very
    specific multi-scale geometry relevant to some models of Active Galacic Nuclei (AGN).

    <b>Background and motivation</b>

    The AGN unification model postulates that the central engine is surrounded by an obscuring
    torus of cold gas and dust which reprocesses the primary radiation. Sudden column density
    changes on timescales of hours to days reveal that this obscurer is clumpy, with individual
    clumps distributed across two orders of magnitude in radius. Recent models (Buchner et al.
    2019, Leftley et al. 2024) capture this geometry with many thousands of small, spherical
    clumps, following distributions constrained by observational variability studies.

    This introduces a severe multi-scale problem for the spatial discretization in SKIRT. The torus
    spans parsec scales while individual clumps are orders of magnitude smaller. The commonly used
    grids cannot simultaneously resolve both. Standard adaptive grids (octrees, k-d trees) would
    require an intractable number of cells, or fail to resolve individual clumps. The current class
    offers an alternative, dedicated spatial grid that follows the exact geometry of the clumps and
    enables efficient radiation transfer simulations.

    <b>%Geometry of the grid</b>

    The spatial domain of this grid has the form of a torus centered on the origin, with
    half-opening angle \f$0<\delta<\pi/2\f$ and inner and outer radii \f$r_\text{min}>0\f$ and
    \f$r_\text{max}>r_\text{min}\f$. The domain is thus bounded by a cone and by two concentric
    spheres. Any medium outside of the domain is ignored.

    The grid cells partitioning the spatial domain consist of a number of spheres with centers and
    radii loaded from an input file, plus a single cell representing the portion of the torus not
    covered by any spheres. The spheres are not allowed to overlap each other, and they must be
    fully contained in the spatial domain without straggling the boundaries. These limitations are
    imposed to avoid ambiguities and to simplify/accelerate the implementation.

    As an important consequence of this setup, medium properties will be (a) uniform within each
    sphere, (b) uniform within the complete uncovered portion of the torus, and (c) nonexistant
    outside of the torus.

    <b>Configuring the medium</b>

    During setup, as for any spatial grid, SKIRT will sample the medium properties within each
    spatial cell from the configured medium distribution. While any type of configuration will
    work, it makes most sense to specify a medium distribution that closely matches the geometry of
    this specialty grid. This can be accomplished by using a ParticleMedium with the
    UniformSmoothingKernel to match the spheres in the grid geometry. The background medium can be
    specified through a TorusGeometry with uniform density or even using a UniformBoxGeometry
    (because the material outside of the spatial domain of the grid is ignored). If the model
    includes multiple media types, different medium components can be combined (for example, one
    ParticleMedium for dust clumps and another one for gas clumps).

    The first four columns of the ParticleMedium import file specify the center position and radius
    of each clump. The same file can be read by this specialty grid class to define the spheres;
    the additional columns in the file are ignored. (If there are multiple ParticleMedium
    instances, the import files can be concatenated to define all spheres).

    The ClumpyTorusSpatialGrid class removes spheres from the imported list that (1) do not fully
    lie inside the torus or (2) overlap any of the previously read spheres. Note that the
    ParticleMedium does NOT do this. This is not a problem for spheres outside the torus because
    any mass in this region will be ignored when sampling the density of the spatial grid cells.
    However, the mass of overlapping spheres (removed from the grid but not from the medium) will
    be included when sampling the grid cell densities, most likely distorting the intended picture.

    */
class ClumpyTorusSpatialGrid : public SpatialGrid
{
    ITEM_CONCRETE(ClumpyTorusSpatialGrid, SpatialGrid, "a specialty spatial grid handling a multi-scale clumpy torus")
        ATTRIBUTE_TYPE_DISPLAYED_IF(ClumpyTorusSpatialGrid, "Level3")

        PROPERTY_DOUBLE(openingAngle, "the half opening angle of the torus")
        ATTRIBUTE_QUANTITY(openingAngle, "posangle")
        ATTRIBUTE_MIN_VALUE(openingAngle, "]0 deg")
        ATTRIBUTE_MAX_VALUE(openingAngle, "90 deg[")

        PROPERTY_DOUBLE(minRadius, "the minimum radius of the torus")
        ATTRIBUTE_QUANTITY(minRadius, "length")
        ATTRIBUTE_MIN_VALUE(minRadius, "]0")

        PROPERTY_DOUBLE(maxRadius, "the maximum radius of the torus")
        ATTRIBUTE_QUANTITY(maxRadius, "length")
        ATTRIBUTE_MIN_VALUE(maxRadius, "]0")

        PROPERTY_STRING(filename, "the name of the file defining the spherical clumps")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This function reads the input file defining the spherical clumps and builds the data
        structures needed for the operation of the grid. */
    void setupSelfBefore() override;

    /** The destructor destroys the custom BVH dat astructure */
    ~ClumpyTorusSpatialGrid();

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the grid, which is 3 for this class. */
    int dimension() const override;

    /** This function returns the number of cells in the grid, which equals the number of spherical
        clumps contained in the torus, plus one for the torus itself. */
    int numCells() const override;

    /** This function returns the bounding box that encloses the grid, i.e. the torus. */
    Box boundingBox() const override;

    /** This function returns the volume of the cell with index \f$m\f$. For a clump, this is the
        spherical volume. For the torus, this is the volume of the torus minus the combined volume
        of all clumps. */
    double volume(int m) const override;

    /** This function returns the diagonal of the cell with index \f$m\f$. For a clump, this is the
        diagonal of the sphere. For the torus, the function use a simplistic estimate equal to
        twice the distance between the outer and inner radii. */
    double diagonal(int m) const override;

    /** This function returns the index \f$m\f$ of the cell that contains the position
        \f${\bf{r}}\f$. */
    int cellIndex(Position bfr) const override;

    /** This function returns the central location of the cell with index \f$m\f$. For a clump,
        this is the center of the sphere. For the torus, the function first tries a position on the
        x-axis halfway between the inner and outer radii. If this position happens to be inside a
        clump, a random position (in the torus but outside any clumps) is returned. */
    Position centralPositionInCell(int m) const override;

    /** This function returns a random location from the cell with index \f$m\f$. For a clump, a
        random position within the sphere is generated through analytical inversion. For the torus,
        first a random position within the torus is generated through analytical inversion, which
        is then rejected iteratively as long as it happens to be inside one of the clumps. */
    Position randomPositionInCell(int m) const override;

    /** This function creates and hands over ownership of a path segment generator (an instance of
        a PathSegmentGenerator subclass) appropriate for this spatial grid type. */
    std::unique_ptr<PathSegmentGenerator> createPathSegmentGenerator() const override;

protected:
    /** This function writes the intersection of the grid structure with the xy plane to the
        specified SpatialGridPlotFile object. */
    void write_xy(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid structure with the xz plane to the
        specified SpatialGridPlotFile object. */
    void write_xz(SpatialGridPlotFile* outfile) const override;

    /** This function writes the intersection of the grid structure with the yz plane to the
        specified SpatialGridPlotFile object. */
    void write_yz(SpatialGridPlotFile* outfile) const override;

    /** This function writes 3D information for all cells in the grid structure to the specified
        SpatialGridPlotFile object. */
    void write_xyz(SpatialGridPlotFile* outfile) const override;

    //======================== Data Members ========================

private:
    // data type defining a single clump
    class Clump
    {
        double _x, _y, _z, _r;

    public:
        Clump(double x, double y, double z, double r) : _x{x}, _y{y}, _z{z}, _r{r} {}
        Position center() const { return Position(_x, _y, _z); }
        double radius() const { return _r; }
        Box bounds() const { return Box(_x - _r, _y - _r, _z - _r, _x + _r, _y + _r, _z + _r); }
    };

    // array defining the clumps
    int _numClumps{0};      // the number of clumps AND the index of the cell representing the torus
    vector<Clump> _clumps;  // index on m, assuming m < _numClumps

    // the custom bounding volume hierarchy that allows efficient querying for our purposes
    class BVH;
    BVH* _bvh{nullptr};

    // allow our path segment generator to access our private data members
    class MySegmentGenerator;
    friend class MySegmentGenerator;
};

////////////////////////////////////////////////////////////////////

#endif
