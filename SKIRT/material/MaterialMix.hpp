/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MATERIALMIX_HPP
#define MATERIALMIX_HPP

#include "Array.hpp"
#include "Direction.hpp"
#include "MaterialMixProvider.hpp"
#include "UpdateStatus.hpp"
class Configuration;
class HistoryScope;
class MaterialState;
class MediumState;
class StateVariableList;
class PhotonPacket;
class Random;
class StokesVector;
class WavelengthGrid;
class DisjointWavelengthGrid;

////////////////////////////////////////////////////////////////////

/** MaterialMix is the abstract base class for all classes representing the concrete material
    properties of a specific transfer medium. The MaterialMix class hierarchy allows fundamentally
    different material types (e.g. dust, electrons, and hydrogen-dominated gas) to be implemented
    as part of a single framework.

    Instances of MaterialMix subclasses are immutable after setup has been completed, so the same
    instance can be reused in multiple contexts.

    <b>Material properties</b>

    The medium state maintained by a simulation for each cell and medium component includes a
    pointer to a MaterialMix instance defining the properties of the material, and a number density
    value defining the amount of material present in the cell per unit of volume. The kind of
    physical entity being counted by the number density and the conversion from number density to
    mass density depend on the type of material, as indicated in the table below.

    Material type | Entity counted | Mass conversion
    --------------|----------------|-----------------------------
    Dust          | hydrogen atom  | dust mass per hydrogen atom
    Electrons     | electron       | electron mass
    Gas           | hydrogen atom  | gas mass per hydrogen atom

    The following table lists some relevant physical quantities including cell properties that may
    be traced by a simulation, material properties defined by material mixes, and properties that
    can be derived from these.

    <TABLE>
    <TR><TD><B>Symbol</B></TD>  <TD><B>Units</B></TD>  <TD><B>Description</B></TD></TR>
    <TR><TD>\f$\Delta s\f$</TD>  <TD>\f$m\f$</TD>  <TD>Distance along a path</TD></TR>
    <TR><TD>\f$V\f$</TD>  <TD>\f$\text{m}^3\f$</TD>  <TD>Volume</TD></TR>
    <TR><TD>\f$v\f$</TD>  <TD>\f$\text{m}\,\text{s}^{-1}\f$</TD>  <TD>Bulk velocity</TD></TR>
    <TR><TD>\f$\bf{B}\f$</TD>  <TD>\f$\text{T}\f$</TD>  <TD>Magnetic field vector</TD></TR>
    <TR><TD>\f$T\f$</TD>  <TD>\f$\text{K}\f$</TD>  <TD>Temperature</TD></TR>
    <TR><TD>\f$n\f$</TD>  <TD>\f$\#\,\text{m}^{-3}\f$</TD>  <TD>Number density (of entities)</TD></TR>
    <TR><TD>\f$\mu\f$</TD>  <TD>\f$\text{kg}\,\#^{-1}\f$</TD>  <TD>Mass per entity</TD></TR>
    <TR><TD>\f$\varsigma\f$</TD>  <TD>\f$\text{m}^2\,\#^{-1}\f$</TD>  <TD>Cross section per entity</TD></TR>
    <TR><TD>\f$\mathcal{N}=n\Delta s\f$</TD> <TD>\f$\#\,\text{m}^{-2}\f$</TD>  <TD>Number column density</TD></TR>
    <TR><TD>\f$N=nV\f$</TD>  <TD>\f$\#\f$</TD>  <TD>Number (of entities)</TD></TR>
    <TR><TD>\f$\rho=n\mu\f$</TD>  <TD>\f$\text{kg}\,\text{m}^{-3}\f$</TD>  <TD>Mass density</TD></TR>
    <TR><TD>\f$\Sigma=n\mu\Delta s\f$</TD> <TD>\f$\text{kg}\,\text{m}^{-2}\f$</TD>  <TD>Mass column density</TD></TR>
    <TR><TD>\f$M=n\mu V\f$</TD>  <TD>\f$\text{kg}\f$</TD>  <TD>Mass</TD></TR>
    <TR><TD>\f$\kappa=\varsigma/\mu\f$</TD>  <TD>\f$\text{m}^2\,\text{kg}^{-1}\f$</TD>  <TD>Mass coefficient</TD></TR>
    <TR><TD>\f$k=n\varsigma\f$</TD>  <TD>\f$\text{m}^{-1}\f$</TD>  <TD>Opacity</TD></TR>
    <TR><TD>\f$\tau=n\varsigma\Delta s\f$</TD>  <TD>\f$1\f$</TD>  <TD>Optical depth</TD></TR>
    </TABLE>

    <b>Capabilities functions</b>

    The abstract interface defined by the MaterialMix base class allows handling essentially all
    aspects of the simulation with respect to a material of a given type, including absorption,
    scattering, and secondary emission of photon packets. The base class does not provide any
    functionality, so everything must be implemented in subclasses.

    Because not all physical processes are relevant for (or supported by) all material types, the
    MaterialMix interface includes a set of functions to define the capabilities of a given
    concrete material mix. On a basic level, the materialType() function returns the overall
    material category (dust, gas, or electrons). In addition, there are a number of Boolean
    functions that indicate whether a certain physical process is supported.

    This approach allows fine-grained run-time discovery of capabilities. The functions can be
    used, for example, during setup to ensure that the configuration is valid (e.g., all material
    mixes support stochastic heating when enabled in the configuration), to enable or disable
    optimizations (e.g., when calculating optical depth), and to enable probing of the appropriate
    information (e.g., producing separate density cuts for dust, gas, and electrons).

    The capabilities() function, declared in the MaterialMixProvider base class, gathers the
    answers of all these functions in a single structure. It is implemented once in this class,
    for all material mixes, so that a subclass overrides the individual functions only. The code
    that establishes the configuration of a simulation asks each medium component for this
    structure (see Medium::capabilities()), so that it does not need a material mix for a medium
    component that is configured with a material mix family.

    Most MaterialMix subclasses implement no public functions outside of those defined in this base
    class. However, sometimes specific features require external access to additional information
    offered by a material mix. For example, MultiGrainDustMix subclasses offer access to individual
    grain populations for use in probes or dust destruction recipes. In those contexts, it is
    acceptable to apply a dynamic cast to discover whether a given material mix offers the relevant
    extended interface.

    <b>%Medium state setup functions</b>

    The MaterialMix class hierarchy offers a number of functions that advertise the required medium
    state variables and assist with initializing their values during setup. For example, the
    specificStateVariableInfo() function returns a list of medium state variable descriptors
    specifying the specific state variables used by the material mix. This allows the medium system
    to allocate storage for the appropriate set of state variables. A material mix with custom
    state variables assembles this list during setup through a StateVariableList, which hands out a
    handle for each variable (see the StateScalar and StateArray classes); the material mix then
    uses these handles to access the variables in the medium state.

    All common state variables and the number density (part of the specific state) are initialized
    by the medium system. Additional specific state variables must be initialized in the
    initializeSpecificState() material mix function, which is invoked by the medium system for each
    spatial cell just after the common state variables and the number density have been
    initialized. If the material mix is configured as part of an imported medium component, extra
    data fields imported from the snapshot based on the information returned by the parameterInfo()
    function are passed to this function. A material mix that needs the temperature or the
    metallicity of the material requests these quantities through the parameterInfo() function,
    like any other quantity, and stores them in custom state variables with the standard names \c
    T and \c Z (see the StateVariable class).

    <b>%Medium state update functions</b>

    As described in the documentation for the MonteCarloSimulation class, the simulation execution
    flow supports a dynamic \em medium \em state (DMS). Specifically, a material mix can provide an
    algorithm for updating its specific state at the end of a primary or secondary emission segment
    based on the radiation field calculated during that segment. The hasDynamicMediumState()
    function indicates if and when the state for this material mix should be updated, and the
    updateSpecificState() function actually performs such an update.

    A material mix can also perform the update in four stages, as indicated by the
    hasStagedStateUpdate() function. Stage 1, the prepareSpecificStateUpdate() function, is called
    for each spatial cell in parallel, and may store in the specific state of the cell the
    information needed by the next stage. After the medium state has been synchronized between
    processes, stage 2, the beginSpecificStateUpdate() function, is called once on each process,
    serially, with access to the state of all cells. Stage 3, the updateSpecificState() function,
    is then called for each cell in parallel, as for any material mix with a dynamic medium state,
    and after the next synchronization, stage 4, the endSpecificStateUpdate() function, is called
    once on each process, serially. Because the state has been synchronized, all processes see the
    same state in stages 2 and 4, and they arrive at the same results without communicating.
    Stage 2 may store data for use by stage 3 in \c mutable data members of the material mix, and
    stage 4 releases them. This is safe because stages 2 and 4 are serial, and stage 3 only reads
    the data. The lifetime of such data is limited to a single update; anything that must survive
    an update belongs in the medium state. Stage 2 may itself perform a parallel calculation
    through the ParallelFactory, and communicate the results between processes through the
    ProcessManager. The TrivialStagedGasMix class illustrates the protocol.

    <b>Low-level material properties functions</b>

    The MaterialMix class hierarchy offers functions for retrieving some basic material properties
    as a function of wavelength, including the absorption cross section, the scattering cross
    section, and the scattering asymmetry parameter. These functions return \em default property
    values, assuming fixed, predefined values for any quantities other than wavelength (e.g., a
    default temperature, no polarization, no kinematics).

    In principle, the values returned by these low-level functions may be used only during setup
    and for probing. However, some portions of the photon life cycle code might be optimized to use
    these functions directly in cases where the optical properties are known to depend solely on
    the photon packet’s wavelength.

    <b>High-level functions for photon life cycle</b>

    Most importantly, the MaterialMix class hierarchy offers a set of functions that help implement
    the photon life cycle on a high, generic level. These functions receive at least two arguments:
    an object representing the medium state for a spatial cell and for a medium component
    configured with the receiving material mix, and an incoming photon packet. Extra arguments may
    override information that is also available as part of the state or photon packet, or they may
    simply provide additional information.

    For example, the opacityAbs() and opacitySca() functions return the absorption and scattering
    opacity \f$k=n\varsigma\f$. They are given a wavelength that overrides the photon packet
    wavelength. Providing a photon packet is in fact optional so that these functions can be used
    in situations where there is no photon packet involved, such as when calculating the luminosity
    absorbed by the dust in a cell.

    The propagate() function adjusts the photon packet for any effects caused by propagation over a
    given distance through the cell. This may include, for example, changes to the polarization
    state caused by dichroism. The function also returns the total (possibly dichroic) optical
    depth for the photon packet intensity over the given distance.

    The performScattering() function handles a complete random-walk scattering interaction with the
    medium component of the receiving material mix, including the effects of bulk velocity,
    polarization, and so forth. The peelOffScattering() function similarly calculates the
    contribution to a scattering peel-off event for this material, given the instrument reference
    frame and the relative weight of this medium component.

    <b>Functions for secondary continuum emission</b>

    The emissivity() function returns the emissivity spectrum of a material mix when it would be
    embedded in a given radiation field, assuming default values for any specific state variables
    and assuming isotropic, unpolarized emission. For regular dust mixes that don't vary spatially
    and that represent spherical grains, the emissivity tells the full story. For other material
    mixes, it provides just a default.

    The emissionSpectrum() function returns the emission spectrum in the spatial cell and medium
    component represented by the specified material state and the receiving material mix when it
    would be embedded in the specified radiation field. The returned spectrum takes into account
    the values of any relevant specific or common state variables, including the number density of
    the dust in the specified cell.

    The remaining functions in this section provide additional optical properties intended for use
    with the SpheroidalPolarization mode.

    <b>Functions for secondary line emission</b>

    The lineEmissionSpectrum() function returns the luminosities that will be emitted at the line
    centers defined by the lineEmissionCenters() function for a material mix with a given material
    state and embedded in a given radiation field. The caller handling secondary emission will add
    a Doppler shift corresponding to random thermal motion with particle masses given by the
    lineEmissionCenters() function. The intrinsic line widths are ignored because they are usually
    much smaller than the thermal dispersion.

    <b>Indicative temperature function</b>

    The indicativeTemperature() function returns an indicative temperature with an interpretation
    depending on the material type. For dust mixes it returns the averaged equilibrium temperature
    of the grain population given the specified radiation field and assuming local thermal
    equilibrium conditions. Other materials may return a temperature determined based on the
    radiation field, the specific state, a default value, or zero if none of the above apply.
*/
class MaterialMix : public MaterialMixProvider
{
    ITEM_ABSTRACT(MaterialMix, MaterialMixProvider, "a material mix")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function caches the simulation's random generator for use by subclasses. */
    void setupSelfBefore() override;

    //======== Material type =======

public:
    /** This function returns the fundamental material type represented by this material mix. See
        the documentation of the MaterialMix class for more information. The MaterialType
        enumeration is defined in the MaterialMixProvider base class. */
    virtual MaterialType materialType() const = 0;

    /** This convenience function returns true if the fundamental material type represented by this
        material mix is Dust, and false otherwise. */
    bool isDust() const { return materialType() == MaterialType::Dust; }

    /** This convenience function returns true if the fundamental material type represented by this
        material mix is Electrons, and false otherwise. */
    bool isElectrons() const { return materialType() == MaterialType::Electrons; }

    /** This convenience function returns true if the fundamental material type represented by this
        material mix is Gas, and false otherwise. */
    bool isGas() const { return materialType() == MaterialType::Gas; }

    //======== Capabilities =======

public:
    /** This function returns true if this material mix supports polarization during scattering
        events, and false otherwise. The default implementation in this base class returns false.
        */
    virtual bool hasPolarizedScattering() const;

    /** This function returns true if the absorption of radiation for this material mix is dichroic
        (i.e. the absorption cross section depends on the polarization state of incoming photon and
        the polarization state is adjusted during absorption), and false otherwise. If
        hasPolarizedAbsorption() returns true, hasPolarizedScattering() must return true as well.
        The default implementation in this base class returns false. */
    virtual bool hasPolarizedAbsorption() const;

    /** This function returns true if the secondary emission for this material mix is or may be
        polarized and anisotropic, and false otherwise. If hasPolarizedEmission() returns true,
        hasPolarizedScattering() must return true as well. The default implementation in this base
        class returns false. */
    virtual bool hasPolarizedEmission() const;

    /** This function returns true if scattering for this material mix is resonant (such as for
        Lyman-alpha), and false otherwise. The default implementation in this base class returns
        false. */
    virtual bool hasResonantScattering() const;

    /** This function returns true if the extinction cross section (the sum of the absorption and
        scattering cross section) for this material mix can be negative, and false otherwise. The
        default implementation in this base class returns false.

        For all material mixes, the scattering cross section must always be positive or zero for
        all wavelengths and material properties. The absorption cross section at a given wavelength
        can be negative if the material exhibits stimulated emission at that wavelength. As long as
        the magnitude of the absorption cross section is guaranteed to be smaller than the
        scattering cross section, the extinction cross section (the sum of both cross sections)
        always remains positive, and this function can safely return false. As soon as the
        magnitude of the negative absorption cross section can be larger than the scattering cross
        section for some wavelengths and material properties, this function must return true. This
        allows the photon cycle machinery to properly handle negative extinction cross sections
        and the corresponding negative optical depths. */
    virtual bool hasNegativeExtinction() const;

    /** This function returns true if this material mix represents dust and supports stochastic
        heating of dust grains for the calculation of secondary emission, and false otherwise. The
        default implementation in this base class returns false. */
    virtual bool hasStochasticDustEmission() const;

    /** This function returns true if the cross sections returned by this material mix may depend
        on the values of specific state variables other than the number density, and false
        otherwise. The default implementation in this base class returns false. */
    virtual bool hasExtraSpecificState() const;

    /** This function returns true if a scattering interaction for this material mix may adjust the
        wavelength of the interacting photon packet, and false otherwise. The default
        implementation in this base class returns false. */
    virtual bool hasScatteringDispersion() const;

    /** This function returns true if some of the outgoing photon packets of a scattering
        interaction for this material mix may emulate secondary emission, and false otherwise.
        These special photon packets will be recorded by instruments as if they originated during
        secondary emission. This feature is used, for example, to implement fluorescence as a
        scattering interaction. The default implementation in this base class returns false. */
    virtual bool scatteringEmulatesSecondaryEmission() const;

    /** This function returns an enumeration indicating whether the material mix offers an
        algorithm to update its specific medium state, and if so, when the updateSpecificState()
        function should be called. Refer to the documentation of the MonteCarloSimulation class for
        more information on dynamic \em medium \em state (DMS) support. The DynamicStateType
        enumeration is defined in the MaterialMixProvider base class.

        The default implementation in this base class returns DynamicStateType::None. */
    virtual DynamicStateType hasDynamicMediumState() const;

    /** This function returns true if this material mix updates its dynamic medium state in four
        stages, as described in the class header, and false if it updates the state of each cell
        through the updateSpecificState() function only. The function is relevant only if the
        hasDynamicMediumState() function returns anything other than \c None. The default
        implementation in this base class returns false. */
    virtual bool hasStagedStateUpdate() const;

    /** This function returns true if this material mix supports secondary continuum emission, i.e.
        with an emission spectrum discretized over a wavelength grid with adjacent bins. The
        default implementation in this base class returns false. */
    virtual bool hasContinuumEmission() const;

    /** This function returns true if this material mix supports secondary line emission, i.e.
        emission at discrete wavelengths. The default implementation in this base class returns
        false. */
    virtual bool hasLineEmission() const;

    /** This function returns the capabilities of this material mix, gathered from the material
        type and capability functions listed above. Because these functions are virtual, a subclass
        overrides them rather than this function. */
    Capabilities capabilities() const final;

    //======== Medium state setup =======

public:
    /** This function returns the number and type of import parameters required by this particular
        material mix as a list of SnapshotParameter objects. Each of these objects specifies unit
        information and a human-readable descripton for the parameter. The default implementation
        in this base class returns an empty list.

        A material mix that needs the temperature or the metallicity of the material requests these
        quantities through this function, using the SnapshotParameter::temperature() and
        SnapshotParameter::metallicity() descriptors, at the position where it wants to receive
        them. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns false because a material mix does not vary with position. */
    bool hasVariableMix() const override;

    /** This function returns a pointer to this material mix, regardless of the specified values of
        the imported quantities. */
    const MaterialMix* mix(const Array& parameters) override;

    /** This function returns a pointer to this material mix. */
    const MaterialMix* mix() override;

    /** This function returns the list of StateVariable objects describing the specific state
        variables used by the receiving material mix. This allows the MediumSystem class to
        allocate storage for the appropriate set of state variables, and it allows probing the
        relevant medium state variables for output. See the StateVariable and StateVariableList
        classes for more info.

        Common state variables are \em not listed; their presence is derived from other aspects of
        the configured medium components. On the other hand, \em all specific state variables used
        by the material mix are listed, starting with the number density, which is always present.
        A material mix with custom state variables assembles the list during its setup, adding the
        custom variables in order, and keeps it as a data member. The function may therefore be
        called only after the material mix has been set up, and the returned reference remains
        valid as long as the material mix exists.

        The default implementation in this base class returns a list that holds just the number
        density. */
    virtual const StateVariableList& specificStateVariableInfo() const;

    /** This function initializes any specific state variables requested by this material mix
        through the specificStateVariableInfo() function except for the number density. The
        function is invoked by the medium system for each spatial cell after the common state
        variables and the number density have been initialized to their proper values, and the
        specific state variables have been initialized to a default value of zero.

        If the material mix is configured as part of an imported medium component, the quantities
        imported from the snapshot as requested by the parameterInfo() function are passed to this
        function, in the same order. If the material mix is configured in a geometric medium
        component, the array is empty, and the material mix uses default values for these
        quantities, usually provided by its own properties.

        The default implementation in this base class does nothing. */
    virtual void initializeSpecificState(MaterialState* state, const Array& params) const;

    //======== Medium state updates =======

    /** If this material mix has a dynamic medium state, i.e. if the hasDynamicMediumState()
        function returns anything other than \c None, this function is invoked for each spatial
        cell at the end of each relevant primary or secondary emission segment. Based on the
        specified radiation field, if needed, the function updates any values in the specific
        material state for this cell and medium component that may inform the local emission and/or
        extinction properties of the material. The function returns the update status as described
        for the UpdateStatus class. This information is used to optimize the synchronization of
        changes across multiple processes and to determine whether an iterative update has
        converged, if applicable.

        This function may be called in parallel across multiple threads and processes. It should
        not update any information outside the specific medium state of the given cell and medium
        component. For a material mix with a staged update, this function is stage 3 of the update,
        and it may read the data prepared by the beginSpecificStateUpdate() function. The default
        implementation in this base class throws a fatal error. */
    virtual UpdateStatus updateSpecificState(MaterialState* state, const Array& Jv) const;

    /** If this material mix has a staged dynamic medium state update, i.e. if the
        hasStagedStateUpdate() function returns true, this function is invoked for each spatial cell
        as stage 1 of each update, before the updateSpecificState() function. Based on the
        specified radiation field, if needed, the function prepares the update of the cell, and may
        store in the specific state of the cell the information that stage 2 needs. The function
        returns the update status as described for the UpdateStatus class; the state of the cells
        that are marked as updated is synchronized between processes before stage 2.

        This function may be called in parallel across multiple threads and processes. It should
        not update any information outside the specific medium state of the given cell and medium
        component. The default implementation in this base class throws a fatal error. */
    virtual UpdateStatus prepareSpecificStateUpdate(MaterialState* state, const Array& Jv) const;

    /** If this material mix has a staged dynamic medium state update, i.e. if the
        hasStagedStateUpdate() function returns true, this function is invoked as stage 2 of each
        update, after the prepareSpecificStateUpdate() function has been called for all spatial
        cells and the medium state has been synchronized between processes. The function is called
        once on each process, serially, with read access to the complete medium state; the \em h
        argument specifies the index of the medium component configured with this material mix.
        The function may store data for use by the updateSpecificState() function in \c mutable
        data members, as described in the class header. For a medium component with a material
        mix family, the function is invoked on the material mix of the first spatial cell. The
        default implementation in this base class does nothing. */
    virtual void beginSpecificStateUpdate(const MediumState& state, int h) const;

    /** If this material mix has a staged dynamic medium state update, i.e. if the
        hasStagedStateUpdate() function returns true, this function is invoked as stage 4 of each
        update, after the updateSpecificState() function has been called for all spatial cells and
        the medium state has been synchronized between processes, and before the
        isSpecificStateConverged() function. The function is called once on each process, serially,
        with read access to the complete medium state; the \em h argument specifies the index of
        the medium component configured with this material mix. The function should release any
        data stored by the beginSpecificStateUpdate() function. For a medium component with a
        material mix family, the function is invoked on the material mix of the first spatial cell.
        The default implementation in this base class does nothing. */
    virtual void endSpecificStateUpdate(const MediumState& state, int h) const;

    /** If this material mix has a dynamic medium state, i.e. if the hasDynamicMediumState()
        function returns anything other than \c None, this function is invoked (once) \em after
        updateSpecificState() has been called for all spatial cells. The \em numCells, \em
        numUpdated and \em numNotConverged arguments specify respectively the number of spatial
        cells in the simulation, the number of cells updated during this update cycle, and the
        number of updated cells that have not yet converged.

        The \em history argument gives access to the series of the IterationHistory kept on behalf
        of the medium component configured with this material mix. Through this scope, the
        function can retrieve the aggregate series declared for the state variables that the mix
        marked for aggregation (see StateVariable::aggregated()), which hold the volume integral
        of the variable after the current update and after the preceding iterations, and it can
        declare and set further scalar series of its own, such as the fraction of converged cells,
        with ids that it chooses itself. A mix keeps no history in its own data members.

        Based on this information and any relevant user configuration options, the function returns
        true if the medium state is considered to be converged and false if not. The default
        implementation in this base class always returns true. */
    virtual bool isSpecificStateConverged(int numCells, int numUpdated, int numNotConverged,
                                          const HistoryScope& history) const;

    //======== Low-level material properties =======

public:
    /** This function returns the mass per entity \f$\mu\f$ for this material. The table below
        indicates the precise meaning of this number depending on the type of material being
        represented.

        Material type | Interpretation of mass() return value
        --------------|---------------------------------------
        Dust          | dust mass per hydrogen atom
        Electrons     | electron mass
        Gas           | gas mass per hydrogen atom
        */
    virtual double mass() const = 0;

    /** This function returns the default absorption cross section per entity
        \f$\varsigma^{\text{abs}}_{\lambda}\f$ at wavelength \f$\lambda\f$. */
    virtual double sectionAbs(double lambda) const = 0;

    /** This function returns the default scattering cross section per entity
        \f$\varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. */
    virtual double sectionSca(double lambda) const = 0;

    /** This function returns the default extinction cross section per entity
        \f$\varsigma^{\text{ext}}_{\lambda} = \varsigma^{\text{abs}}_{\lambda} +
        \varsigma^{\text{sca}}_{\lambda}\f$ at wavelength \f$\lambda\f$. */
    virtual double sectionExt(double lambda) const = 0;

    /** This function returns the default scattering asymmetry parameter \f$g_\lambda =
        \left<\cos\theta\right>\f$ at wavelength \f$\lambda\f$. This value serves as a parameter
        for the Henyey-Greenstein phase function. The default implementation in this base class
        returns zero, indicating isotropic scattering. */
    virtual double asymmpar(double lambda) const;

    //======== High-level photon life cycle =======

public:
    /** This function returns the absorption opacity \f$k^\text{abs}=n\varsigma^\text{abs}\f$ for
        the given wavelength, material state, and photon properties (optional; may be nullptr). */
    virtual double opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* pp) const = 0;

    /** This function returns the scattering opacity \f$k^\text{sca}=n\varsigma^\text{sca}\f$ for
        the given wavelength, material state, and photon properties (optional; may be nullptr). */
    virtual double opacitySca(double lambda, const MaterialState* state, const PhotonPacket* pp) const = 0;

    /** This function returns the extinction opacity \f$k^\text{ext}=k^\text{abs}+k^\text{sca}\f$
        for the given wavelength, material state, and photon properties (optional; may be nullptr).
        */
    virtual double opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const = 0;

    /** This function calculates the contribution of the medium component associated with this
        material mix to the peel-off photon luminosity, polarization state, and wavelength shift
        for the given wavelength, geometry, material state, and photon properties. The
        contributions to the Stokes vector components are stored in the \em I, \em Q, \em U, \em V
        arguments, which are guaranteed to be initialized to zero by the caller. If there is a
        wavelength shift, the new wavelength value replaces the incoming value of the \em lambda
        argument. The function returns true if the scattering interaction emulates secondary
        emission (e.g., to implement fluorescence as scattering), and false otherwise (i.e. in the
        vast majority of cases).

        Since we force the peel-off photon packet to be scattered from the direction \f${\bf{k}}\f$
        into the direction \f${\bf{k}}_{\text{obs}}\f$, the corresponding biasing factor is given
        by the probability that a photon packet would be scattered into the direction
        \f${\bf{k}}_{\text{obs}}\f$ if its original propagation direction was \f${\bf{k}}\f$. For a
        given medium component, this biasing factor is equal to the value of the scattering phase
        function \f$\Phi({\bf{k}},{\bf{k}}_{\text{obs}})\f$ for that medium component.

        The default implementation in this base class throws a fatal error. */
    virtual bool peeloffScattering(double& I, double& Q, double& U, double& V, double& lambda, Direction bfkobs,
                                   Direction bfky, const MaterialState* state, const PhotonPacket* pp) const;

    /** This function performs a scattering event on the specified photon packet in the spatial
        cell and medium component represented by the specified material state and the receiving
        material mix. Most of the properties of the photon packet remain unaltered, including the
        position and the luminosity. The properties that change include the number of scattering
        events experienced by the photon packet, which is increased by one, the propagation
        direction, which is generated randomly, the wavelength, which is properly Doppler-shifted
        for the bulk velocity of the medium, and the polarization state, which may be affected by
        the scattering process.

        The calculation takes all physical processes into account, including the bulk velocity and
        Hubble expansion velocity in the cell, any relevant material state variables such as the
        temperature of a gas medium, and any relevant properties of the incoming photon packet such
        as the polarization state. The first argument specifies the perceived wavelength of the
        photon packet at the scattering location so that this value does not need to be
        recalculated within the function.

        The default implementation in this base class throws a fatal error. */
    virtual void performScattering(double lambda, const MaterialState* state, PhotonPacket* pp) const;

    //======== Secondary continuum emission =======

public:
    /** If this material mix supports secondary continuum emission, this function returns the
        wavelength grid on which this emission is discretized. For more information, refer to the
        documentation of this function for each material type. The default implementation in this
        base class throws a fatal error. */
    virtual DisjointWavelengthGrid* emissionWavelengthGrid() const;

    /** This function returns the continuum emissivity spectrum \f$\varepsilon_{\ell'}\f$ (radiated
        power per unit of solid angle and per material entity) of the material mix when it would be
        embedded in the radiation field specified by the mean intensities \f$(J_\lambda)_\ell\f$,
        assuming default values for any specific state variables and assuming isotropic,
        unpolarized emission. For regular dust mixes that don't vary spatially and that represent
        spherical grains, the emissivity tells the full story. For other material mixes, it
        provides just a default or it is not implemented.

        The input radiation field must be discretized on the simulation's radiation field
        wavelength grid as returned by the Configuration::radiationFieldWLG() function. The output
        emissivity spectrum is discretized on the wavelength grid returned by the
        emissionWavelengthGrid() function. For more information, refer to the documentation of this
        function for each material type. The default implementation in this base class throws a
        fatal error. */
    virtual Array emissivity(const Array& Jv) const;

    /** This function returns the continuum emission spectrum (radiated power per unit of solid
        angle) in the spatial cell and medium component represented by the specified material state
        and the receiving material mix when it would be embedded in the specified radiation field.
        The returned spectrum takes into account the values of any specific and common state
        variables, including the number density of the material in the specified cell. As a result,
        the spectra returned for multiple media components of the same material type can be
        aggregated through simple summation. However, the caller is responsible for final
        normalization after such aggregation has taken place.

        The input radiation field must be discretized on the simulation's radiation field
        wavelength grid as returned by the Configuration::radiationFieldWLG() function. The output
        emissivity spectrum is discretized on the wavelength grid returned by the
        emissionWavelengthGrid() function. For more information, refer to the documentation of this
        function for each material type. The default implementation in this base class throws a
        fatal error. */
    virtual Array emissionSpectrum(const MaterialState* state, const Array& Jv) const;

    /** This function is intended for use with the SpheroidalPolarization mode. It returns the grid
        used for discretizing quantities that are a function of the scattering/emission angle
        \f$\theta\f$. The same grid is returned by all material mixes that have
        SpheroidalPolarization mode. The default implementation in this base class throws a fatal
        error. */
    virtual const Array& thetaGrid() const;

    /** This function is intended for use with the SpheroidalPolarization mode. It returns the
        absorption cross sections per entity \f$\varsigma ^{\text{abs}} _{\lambda} (\theta)\f$ at
        wavelength \f$\lambda\f$ as a function of the emission angle \f$\theta\f$, discretized on
        the grid returned by the thetaGrid() function. The default implementation in this base
        class throws a fatal error. */
    virtual const Array& sectionsAbs(double lambda) const;

    /** This function is intended for use with the SpheroidalPolarization mode. It returns the
        linear polarization absorption cross sections per entity \f$\varsigma ^{\text{abspol}}
        _{\lambda} (\theta)\f$ at wavelength \f$\lambda\f$ as a function of the emission angle
        \f$\theta\f$, discretized on the grid returned by the thetaGrid() function. The default
        implementation in this base class throws a fatal error. */
    virtual const Array& sectionsAbspol(double lambda) const;

    //======== Secondary line emission =======

public:
    /** If this material mix supports secondary line emission, this function returns a list of the
        line centers at which this emission occurs. The default implementation in this base class
        throws a fatal error. */
    virtual Array lineEmissionCenters() const;

    /** If this material mix supports secondary line emission, this function returns a list of the
        particle masses for each of the lines returned by the lineEmissionCenters() function. These
        masses can be used by a caller in conjunction with the temperature stored in the medium
        state to determine the thermal velocity distribution of emitting particles and assign a
        corresponding random Doppler shift to the wavelength of emitted photon packets. The default
        implementation in this base class throws a fatal error. */
    virtual Array lineEmissionMasses() const;

    /** This function returns the line emission spectrum (radiated power per unit of solid angle)
        in the spatial cell and medium component represented by the specified material state and
        the receiving material mix when it would be embedded in the specified radiation field. The
        returned spectrum takes into account the values of any specific and common state variables,
        including the number density of the material in the specified cell. The returned values
        correspond to each of the lines returned by the lineEmissionCenters() function. The default
        implementation in this base class throws a fatal error. */
    virtual Array lineEmissionSpectrum(const MaterialState* state, const Array& Jv) const;

    //============== Indicative temperature =============

public:
    /** This function returns an indicative temperature for the material represented by the
        specified material state and the receiving material mix, assuming an embedding radiation
        field specified by the mean intensities \f$(J_\lambda)_\ell\f$, if available.

        If the simulation tracks the radiation field, the specified \em Jv array is discretized on
        the simulation's radiation field wavelength grid as returned by the
        Configuration::radiationFieldWLG() function. If the simulation does not track the radiation
        field, the \em Jv array passed to this function is empty.

        The interpretation of the indicative temperature depends heavily on the material type. For
        dust mixes, the function returns the averaged equilibrium temperature of the grain
        population given the specified radiation field and assuming local thermal equilibrium
        conditions. Other materials may return a temperature determined based on the radiation
        field and/or the material state, a default value, or zero if none of the above apply. Refer
        to the description for this function in the various subclasses for more information. The
        default implementation in this base class throws a fatal error. */
    virtual double indicativeTemperature(const MaterialState* state, const Array& Jv) const;

    //======================== Other Functions =======================

protected:
    /** This function returns the simulation's random generator as a service to subclasses. */
    Random* random() const { return _random; }

    /** This function returns the simulation's configuration object as a service to subclasses. */
    Configuration* config() const { return _config; }

    //======================== Data Members ========================

private:
    // data members initialized in setupSelfBefore
    Random* _random{nullptr};
    Configuration* _config{nullptr};
};

////////////////////////////////////////////////////////////////////

#endif
