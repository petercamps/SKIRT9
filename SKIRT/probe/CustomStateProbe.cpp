/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CustomStateProbe.hpp"
#include "ColumnOutFile.hpp"
#include "Configuration.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "ProbeFormBridge.hpp"
#include "StateVariableList.hpp"
#include "StringUtils.hpp"
#include "Units.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // the information on a single custom value, i.e. a scalar custom state variable or an element of an array variable
    struct CustomValue
    {
        string name;         // the name of the variable, if any
        string description;  // the description of the variable, or of the element of an array variable
        string quantity;     // the physical quantity represented by the value
        char format;         // the output format specifier
        int customIndex;     // the custom index of the value
    };

    // returns the information on the custom values of the specified list of state variables, in order of custom index
    vector<CustomValue> listCustomValues(const StateVariableList& variables)
    {
        vector<CustomValue> result;
        for (const StateVariable& variable : variables)
        {
            if (variable.identifier() != StateVariable::Identifier::Custom) continue;
            for (int k = 0; k != variable.count(); ++k)
                result.push_back({variable.name(), variable.elementDescription(k), variable.quantity(),
                                  variable.format(), variable.customIndex() + k});
        }
        return result;
    }

    // returns a list of integers corresponding to the names, indices, and index ranges in the specified string,
    // given the information on the custom values; the name of an array variable selects all its elements
    // - if the string is empty, all indices 0...numIndices are included
    // - otherwise, the string must be a comma-separated list of names, zero-based indices, or index ranges
    // portions of the string that don't conform to the syntax, and unknown names, are ignored
    vector<int> listIndices(const vector<CustomValue>& descriptors, string list)
    {
        int n = descriptors.size();
        vector<int> result;
        list = StringUtils::squeeze(list);
        if (list.empty())
        {
            for (int i = 0; i != n; ++i) result.push_back(i);
        }
        else
        {
            for (string segment : StringUtils::split(list, ","))
            {
                segment = StringUtils::squeeze(segment);
                auto parts = StringUtils::split(segment, "-");
                if (parts.size() == 1 && !StringUtils::isValidInt(parts[0]))
                {
                    for (int i = 0; i != n; ++i)
                        if (descriptors[i].name == segment) result.push_back(i);
                }
                else if (parts.size() == 1)
                {
                    int index = StringUtils::toInt(parts[0]);
                    if (index >= n) continue;
                    result.push_back(index);
                }
                else if (parts.size() == 2)
                {
                    if (!StringUtils::isValidInt(parts[0])) continue;
                    if (!StringUtils::isValidInt(parts[1])) continue;
                    int index1 = StringUtils::toInt(parts[0]);
                    int index2 = StringUtils::toInt(parts[1]);
                    if (index2 >= n) index2 = n - 1;
                    for (int i = index1; i <= index2; ++i) result.push_back(i);
                }
            }
        }
        return result;
    }
}

////////////////////////////////////////////////////////////////////

void CustomStateProbe::probe()
{
    if (find<Configuration>()->hasMedium())
    {
        // locate the medium system and units system
        auto ms = find<MediumSystem>();
        auto units = find<Units>();

        // construct a bridge
        ProbeFormBridge bridge(this, form());

        // loop over the media
        int numMedia = ms->numMedia();
        for (int h = 0; h != numMedia; ++h)
        {
            // get the information on the custom values for this medium
            auto descriptors = listCustomValues(ms->media()[h]->mix()->specificStateVariableInfo());

            // get the indices corresponding to the user property and to the nr of custom variables for this medium
            auto indexv = listIndices(descriptors, variables());

            // if there are custom variables to be probed for this medium
            if (!indexv.empty())
            {
                // construct the "quantity axis" from the zero-based indices
                Array axis = NR::array(indexv);

                // construct a list of corresponding unit conversion factors (for performance reasons)
                Array conv(axis.size());
                int outindex = 0;
                for (int index : indexv)
                {
                    const auto& descriptor = descriptors[index];
                    conv[outindex] = descriptor.quantity.empty() ? 1. : units->out(descriptor.quantity, 1.);
                    outindex++;
                }

                // define the call-back function to add column definitions
                auto addColumnDefinitions = [&indexv, &descriptors, units](ColumnOutFile& outfile) {
                    for (int index : indexv)
                    {
                        const auto& descriptor = descriptors[index];
                        string qty = descriptor.quantity;
                        string unit = qty.empty() ? "1" : units->unit(qty);
                        outfile.addColumn(descriptor.description, unit, descriptor.format);
                    }
                };

                // define the call-back function to retrieve a compound value in output units
                auto valueInCell = [ms, &indexv, &descriptors, h, conv](int m) {
                    Array out(indexv.size());
                    int outindex = 0;
                    for (int index : indexv)
                    {
                        const auto& descriptor = descriptors[index];
                        out[outindex] = conv[outindex] * ms->custom(m, h, descriptor.customIndex);
                        outindex++;
                    }
                    return out;
                };

                // define the call-back function to retrieve a cell weight; the units don't matter
                auto weightInCell = [ms](int m) { return ms->massDensity(m); };

                // get the units of the first variable being probed
                string qty = descriptors[indexv[0]].quantity;
                string unit = qty.empty() ? "1" : units->unit(qty);

                // output the requested variables for this medium
                string sh = std::to_string(h);
                bridge.writeQuantity(sh + "_customstate", unit, "custom state quantities",
                                     "density-weighted custom state quantities", axis, "1", addColumnDefinitions,
                                     valueInCell, weightInCell);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
