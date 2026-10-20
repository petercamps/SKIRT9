/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TrivialStagedGasMix.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "MaterialState.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProcessManager.hpp"
#include <set>

////////////////////////////////////////////////////////////////////

void TrivialStagedGasMix::setupSelfBefore()
{
    TrivialGasMix::setupSelfBefore();

    _node = _variables.add(StateVariable::custom("node index", "", 'd').named("node").categorical());
    _nodeDensity = _variables.add(StateVariable::custom("node number density", "numbervolumedensity").named("n_node"));
    _numNodes = _variables.add(StateVariable::custom("number of nodes", "", 'd').named("numNodes"));
}

////////////////////////////////////////////////////////////////////

MaterialMix::DynamicStateType TrivialStagedGasMix::hasDynamicMediumState() const
{
    return DynamicStateType::Primary;
}

////////////////////////////////////////////////////////////////////

bool TrivialStagedGasMix::hasStagedStateUpdate() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

const StateVariableList& TrivialStagedGasMix::specificStateVariableInfo() const
{
    return _variables;
}

////////////////////////////////////////////////////////////////////

UpdateStatus TrivialStagedGasMix::prepareSpecificStateUpdate(MaterialState* state, const Array& /*Jv*/) const
{
    UpdateStatus status;
    double n = state->numberDensity();
    if (n > 0.)
    {
        _node(state) = std::floor(log10(n) / nodeSpacing());
        status.updateConverged();
    }
    return status;
}

////////////////////////////////////////////////////////////////////

void TrivialStagedGasMix::beginSpecificStateUpdate(const MediumState& state, int h) const
{
    // collect the distinct node indices of the cells with material
    std::set<int> nodes;
    for (int m = 0; m != state.numCells(); ++m)
    {
        const MaterialState mst(state, m, h);
        if (mst.numberDensity() > 0.) nodes.insert(static_cast<int>(_node(&mst)));
    }
    _tableNodes.assign(nodes.begin(), nodes.end());

    // calculate the number density represented by each node, distributing the nodes over threads and processes,
    // and collect the results on all processes
    int numNodes = _tableNodes.size();
    _tableDensities.resize(numNodes);
    find<ParallelFactory>()->parallelDistributed()->call(numNodes, [this](size_t firstIndex, size_t numIndices) {
        for (size_t i = firstIndex; i != firstIndex + numIndices; ++i)
            _tableDensities[i] = pow(10., _tableNodes[i] * nodeSpacing());
    });
    ProcessManager::sumToAll(_tableDensities);
}

////////////////////////////////////////////////////////////////////

UpdateStatus TrivialStagedGasMix::updateSpecificState(MaterialState* state, const Array& /*Jv*/) const
{
    UpdateStatus status;
    if (state->numberDensity() > 0.)
    {
        int node = static_cast<int>(_node(state));
        auto position = std::lower_bound(_tableNodes.begin(), _tableNodes.end(), node);
        if (position == _tableNodes.end() || *position != node)
            throw FATALERROR("The node of a spatial cell is missing from the table built in stage 2");
        _nodeDensity(state) = _tableDensities[position - _tableNodes.begin()];
        _numNodes(state) = _tableNodes.size();
        status.updateConverged();
    }
    return status;
}

////////////////////////////////////////////////////////////////////

void TrivialStagedGasMix::endSpecificStateUpdate(const MediumState& /*state*/, int h) const
{
    find<Log>()->info("  TrivialStagedGasMix of medium component " + std::to_string(h) + " used a table with "
                      + std::to_string(_tableNodes.size()) + " nodes");
    _tableNodes.clear();
    _tableDensities.resize(0);
}

////////////////////////////////////////////////////////////////////
