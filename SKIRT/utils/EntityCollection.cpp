/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "EntityCollection.hpp"
#include "Vec.hpp"

////////////////////////////////////////////////////////////////////

EntityCollection::EntityCollection() {}

////////////////////////////////////////////////////////////////////

void EntityCollection::clear()
{
    _entities.clear();
}

////////////////////////////////////////////////////////////////////

void EntityCollection::add(int m, double w)
{
    if (m >= 0 && w > 0. && std::isfinite(w)) _entities.emplace(m, w);
}

////////////////////////////////////////////////////////////////////

void EntityCollection::addSingle(int m)
{
    _entities.clear();
    if (m >= 0) _entities.emplace(m, 1.);
}

////////////////////////////////////////////////////////////////////

double EntityCollection::accumulate(const std::function<double(int)>& value)
{
    double sumvw = 0.;
    for (const auto& [m, w] : _entities)
    {
        sumvw += value(m) * w;
    }
    return sumvw;
}

////////////////////////////////////////////////////////////////////

std::pair<double, double> EntityCollection::average(const std::function<double(int m)>& value,
                                                    const std::function<double(int m)>& weight)
{
    double sumvw = 0.;
    double sumw = 0.;
    for (const auto& [m, wm] : _entities)
    {
        double v = value(m);
        double w = weight(m) * wm;
        sumvw += v * w;
        sumw += w;
    }
    return std::make_pair(sumvw, sumw);
}

////////////////////////////////////////////////////////////////////

double EntityCollection::averageValue(std::function<double(int)> value, const std::function<double(int)>& weight)
{
    auto numEntities = _entities.size();
    if (numEntities == 0) return 0.;
    if (numEntities == 1) return value(_entities.cbegin()->first);

    auto [sumvw, sumw] = average(value, weight);
    return sumw > 0. ? sumvw / sumw : 0.;
}

////////////////////////////////////////////////////////////////////

std::pair<Vec, double> EntityCollection::average(const std::function<Vec(int m)>& value,
                                                 const std::function<double(int m)>& weight)
{
    Vec sumvw;
    double sumw = 0.;
    for (const auto& [m, wm] : _entities)
    {
        Vec v = value(m);
        double w = weight(m) * wm;
        sumvw += v * w;
        sumw += w;
    }
    return std::make_pair(sumvw, sumw);
}

////////////////////////////////////////////////////////////////////

Vec EntityCollection::averageValue(std::function<Vec(int)> value, const std::function<double(int)>& weight)
{
    auto numEntities = _entities.size();
    if (numEntities == 0) return Vec();
    if (numEntities == 1) return value(_entities.cbegin()->first);

    auto [sumvw, sumw] = average(value, weight);
    return sumw > 0. ? sumvw / sumw : Vec();
}

////////////////////////////////////////////////////////////////////
