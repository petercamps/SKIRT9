/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MEDIUMSTATE_HPP
#define MEDIUMSTATE_HPP

#include "Array.hpp"
#include "Span.hpp"
#include "StateVariable.hpp"
#include "UpdateStatus.hpp"
#include "Vec.hpp"
class StateVariableList;

//////////////////////////////////////////////////////////////////////

/** An instance of the MediumState class contains the complete medium state for a simulation. The
    MediumSystem maintains a single (private) instance of this class and provides access to its
    contents to the rest of the simulation code as needed.

    For each spatial cell in the simulation, the medium state includes a set of \em common state
    variables shared by all medium components and a set of \em specific state variables for each
    individual medium component. The following table lists the variables that may be present in
    each of these sets.

    | Identifier | Symbol | Type | Set | Present
    | -----------|--------|------|-----|--------
    | Volume        | \f$V\f$      | double | Common   | always
    | BulkVelocity  | \f$v\f$      | Vec    | Common   | if needed for any medium component
    | MagneticField | \f$\bf{B}\f$ | Vec    | Common   | if needed for any medium component
    | NumberDensity | \f$n\f$      | double | Specific | always
    | Custom + k    |              | double | Specific | zero or more as requested by the medium component

    The temperature and the metallicity of a medium component, if any, are custom variables with
    the standard names \c T and \c Z (see the StateVariable class).

    The presence or absence of these state variables is communicated to the MediumState class
    during its construction, so that it can dynamically allocate the appropriate amount of storage
    and provide an access mechanism for each state variable. This process is described in more
    detail in the following sections.

    <b>Construction</b>

    Fully initializing a MediumState object requires calling the following functions in the
    correct order (after default-constructing the object):
     - the initConfiguration() function specifies the number of spatial cells and the number of
       medium components in the simulation.
     - the initCommonStateVariables() function specifies the set of common state variables.
     - the initSpecificStateVariables() function must be called once for each medium component, in
       order of component index, specifying the set of specific state variables for that component.
     - the initAllocate() function finalizes construction and actually allocates storage;
       it initializes all variables to a value of zero.
     - the setXXX() functions set any nonzero initial variable values required to reflect the input
       model; this may happen in parallel.
     - the initCommunicate() function communicates the initialized state variable values between
       processes.

    The initCommonStateVariables() function receives a list of StateVariable objects to identify
    the required common state variables, including those that should always be present. The
    initSpecificStateVariables() function receives the StateVariableList of the material mix of the
    medium component, which starts with the number density and assigns consecutive custom indices
    in the range \f$ 0 \le k < K\f$, where K is the total number of custom values of the
    component (an array variable occupies as many indices as it has elements). The custom
    variables of a component are stored contiguously for each spatial cell, so that the customs()
    function can return a view on a range of them.

    <b>Synchronization</b>

    A simulation may update the values of medium state variables, often as part of an iterative
    process to calculate a self-consistent state. These updates may happen in parallel over cells.
    After a set of updates, usually at the end of each step in the iterative process, the
    synchronize() function should be called to broadcast the updates between processes.

    <b>Aggregation</b>

    When a simulation employs an iterative process to calculate a self-consistent medium state, it
    needs a criterion to determine whether the state has sufficiently converged. Such a criterion
    is often based on aggregate information calculated over all cells in the spatial grid, such as
    the volume integral of a density. The volumeIntegrals() function calculates the volume
    integrals of a set of custom variables in a single pass over the cells. The medium system
    records these values in the aggregate series of the IterationHistory, which keeps them for the
    current and previous iterations.

    <b>Storage</b>

    To simplify the storage mechanism, state variables of type \c Vec are split into their three
    components, so that all state variables can considered to be of type \c double. Memory is
    allocated only for state variables that are actually needed. Given the number of spatial cells
    \f$M\f$, the number of medium
    components \f$H\f$, the number of common state variables \f$C\f$, and the number of specific
    state variables \f$S_h\f$ for each medium component \f$h\f$, the number of state variables per
    spatial cell is \f$K = C+\sum_h S_h\f$ and the grand total number of state variables is \f$N =
    M (C+\sum_h S_h)\f$.

    The MediumState class allocates a single one-dimensional data array of this size \f$N\f$ and
    provides a mapping to locate a particular state variable in this array. This is accomplished by
    also storing the offset for each variable within the set of variables per spatial cell. To
    facilitate access, these offsets are stored for all supported variables (by name) and not just
    for the required variables. The offsets for unused variables remain at zero. The array is a
    standard vector, so that capacity for additional cells can be reserved in advance.

    Given the offset \f$O_x\f$ for a particular state variable \f$x\f$, a spatial cell index
    \f$m\f$ and a medium component index \f$h\f$, the index of the variable in the data array can
    be calculated as \f$i=K \times m + O_x\f$. This implies that variables are stored contiguously
    per cell.

    <b>Dynamic grid refinement</b>

    When dynamic grid refinement subdivides a cell, the first child keeps the cell index of its
    parent, and the other children are appended to the cell list. The appendCells() function grows
    the medium state accordingly, copying the state of the parent to each new cell, and the
    resizeCell() function then sets the volume of each child, scaling its extensive variables by
    the ratio of the new to the old volume. As a result, intensive and categorical variables are
    copied to each child, and extensive variables are divided among the children in proportion to
    their volume (see the QuantityKind enumeration). The common state variables and the number
    density are intensive.

    <b>Access to undefined variables</b>

    In general, an attempt to access (read or write) a variable for which storage has not been
    requested results in undefined behavior. However, assuming that storage has been requested for
    the volume (which is always needed anyway), the other variables in the common state can be
    safely read (not written) even if no storage was requested for them. This allows the bulk
    velocity and magnetic field vectors to be retrieved without concern for whether they were
    configured in the input model. */
class MediumState
{
    //============= Construction =============

public:
    /** This function initializes the number of spatial cells and number of medium components. */
    void initConfiguration(int numCells, int numMedia);

    /** This function initializes the set of required common state variables. */
    void initCommonStateVariables(const vector<StateVariable>& variables);

    /** This function initializes the set of required specific state variables for the next medium
        component, given the list of variables of its material mix. The function must be called
        once for each medium component, in order of component index. */
    void initSpecificStateVariables(const StateVariableList& variables);

    /** This function ends the initialization sequence, allocates memory for the state variables,
        and returns the total number \f$N = M (C+\sum_h S_h)\f$ of state variables. All newly
        allocated state variables are guaranteed to have a value of zero. */
    size_t initAllocate();

    /** This function communicates the state variable values between processes after each process
        has initialized the values for a subset of the spatial cells and left the values for the
        other cells at zero. (The function uses the ProcessManager::sumToAll() function, so it
        assumes that the uninitialized variables have a zero value). */
    void initCommunicate();

    //============= Synchronization =============

public:
    /** This function synchronizes the state variables between processes after each process has
        possibly updated some of their values. The provided vector indicates the cells for which
        the calling process has made some changes. The vector must be of length \em numCells as
        passed to the initConfiguration() function, where each element in the vector corresponds to
        the cell with the same index. If the calling process updated the state of a given cell, the
        corresponding vector element must be \em UpdatedConverged or \em UpdatedNotConverged and
        otherwise it must be \em NotUpdated. The function returns a pair of integers specifying the
        total number of updated cells and the total number of not-converged cells, aggregated over
        all processes.

        Only one of the calling processes may have updated the state for any given cell. If two or
        more processes updated the state of the same cell, the result of the synchronization is
        undefined. */
    std::pair<int, int> synchronize(const vector<UpdateStatus>& cellFlags);

    //============= Dynamic grid refinement =============

public:
    /** This function grows the medium state by one cell for each element of \em parentCells. The
        new cells receive consecutive indices following the existing cells, and start with the
        values of all state variables of the indicated parent cell, including its volume. The
        function must be called on all processes after the state has been synchronized. */
    void appendCells(const vector<int>& parentCells);

    /** This function sets the volume of the spatial cell with index \f$m\f$ to the specified
        value, and multiplies the extensive specific state variables of the cell by the ratio of the
        new to the old volume. It is intended for the children of a subdivided cell, which start
        with the volume and the state of their parent (see the appendCells() function). */
    void resizeCell(int m, double volume);

    /** This function reserves memory for the specified number of cells, so that the medium state
        can grow to that number of cells without reallocation. */
    void reserveCells(int numCells);

    /** This function returns the number of bytes held by the medium state for each cell. */
    size_t bytesPerCell() const { return _numVars * sizeof(double); }

    //============= Setting =============

public:
    /** This function sets the volume \f$V\f$ of the spatial cell with index \f$m\f$. */
    void setVolume(int m, double value);

    /** This function sets the aggregate bulk velocity \f${\boldsymbol{v}}\f$ of the medium in the
        spatial cell with index \f$m\f$. */
    void setBulkVelocity(int m, Vec value);

    /** This function sets the magnetic field \f${\boldsymbol{B}}\f$ in the spatial cell with index
        \f$m\f$. */
    void setMagneticField(int m, Vec value);

    /** This function sets the number density of the medium component with index \f$h\f$ in the
        spatial cell with index \f$m\f$. */
    void setNumberDensity(int m, int h, double value);

    /** This function sets the value of the custom variable with index \f$i\f$ of the medium
        component with index \f$h\f$ in the spatial cell with index \f$m\f$. */
    void setCustom(int m, int h, int i, double value);

    //============= Querying =============

public:
    /** This function returns the number of spatial cells. */
    int numCells() const { return _numCells; }

    /** This function returns the volume integral over all spatial cells of each of the specified
        custom variables, identified by medium component index and custom variable index,
        calculated in a single pass over the cells. */
    vector<double> volumeIntegrals(const vector<std::pair<int, int>>& variables) const;

    /** This function returns the volume \f$V\f$ of the spatial cell with index \f$m\f$. */
    double volume(int m) const { return _data[_numVars * m + _off_volu]; }

    /** This function returns the aggregate bulk velocity \f${\boldsymbol{v}}\f$ of the medium in
        the spatial cell with index \f$m\f$, or zero if storage was not requested for this
        variable. */
    Vec bulkVelocity(int m) const
    {
        if (_off_velo)
        {
            size_t i = _numVars * m + _off_velo;
            return Vec(_data[i], _data[i + 1], _data[i + 2]);
        }
        return Vec();
    }

    /** This function returns the magnetic field \f${\boldsymbol{B}}\f$ in the spatial cell with
        index \f$m\f$, or zero if storage was not requested for this variable. */
    Vec magneticField(int m) const
    {
        if (_off_mfld)
        {
            size_t i = _numVars * m + _off_mfld;
            return Vec(_data[i], _data[i + 1], _data[i + 2]);
        }
        return Vec();
    }

    /** This function returns the number density of the medium component with index \f$h\f$ in the
        spatial cell with index \f$m\f$. */
    double numberDensity(int m, int h) const { return _data[_numVars * m + _off_dens[h]]; }

    /** This function returns the value of the custom variable with index \f$i\f$ of the medium
        component with index \f$h\f$ in the spatial cell with index \f$m\f$. */
    double custom(int m, int h, int i) const { return _data[_numVars * m + _off_cust[h] + i]; }

    /** This function returns a read-only view on the \f$n\f$ consecutive custom variables starting
        at index \f$i\f$ of the medium component with index \f$h\f$ in the spatial cell with index
        \f$m\f$. */
    Span<const double> customs(int m, int h, int i, int n) const
    {
        return Span<const double>(_data.data() + _numVars * m + _off_cust[h] + i, n);
    }

    /** This function returns a writable view on the \f$n\f$ consecutive custom variables starting
        at index \f$i\f$ of the medium component with index \f$h\f$ in the spatial cell with index
        \f$m\f$. */
    Span<double> customs(int m, int h, int i, int n)
    {
        return Span<double>(_data.data() + _numVars * m + _off_cust[h] + i, n);
    }

    //======================== Data Members ========================

private:
    // data array containing the medium state variables
    vector<double> _data;

    // overall configuration
    int _numCells{0};
    int _numMedia{0};
    size_t _numVars{0};

    // offsets used for mapping common and specific variables (for each medium component) to indices in the data array
    int _off_volu{0};
    int _off_velo{0};
    int _off_mfld{0};
    vector<int> _off_dens;
    vector<int> _off_cust;
    vector<int> _off_extensive;  // offsets of the extensive specific state variables, for all components

    // indices indicating the next item to be initialized; used only during initialization
    int _nextOffset{0};
    int _nextComponent{0};
};

////////////////////////////////////////////////////////////////

#endif
