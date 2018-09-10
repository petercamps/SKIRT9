/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NameManager.hpp"
#include "FatalError.hpp"
#include "BooleanExpression.hpp"

////////////////////////////////////////////////////////////////////

NameManager::NameManager()
{
    void clearAll();
}

////////////////////////////////////////////////////////////////////

void NameManager::clearAll()
{
    _globalSet.clear();
    _globalSet.insert("True");

    while (!_localSetStack.empty()) _localSetStack.pop();
    _localSetStack.emplace(std::initializer_list<string>{{"true"}});
}

////////////////////////////////////////////////////////////////////

void NameManager::pushLocal()
{
    _localSetStack.emplace(std::initializer_list<string>{{"true"}});
}

////////////////////////////////////////////////////////////////////

void NameManager::popLocal()
{
    _localSetStack.pop();
}

////////////////////////////////////////////////////////////////////

namespace
{
    bool isLowercase(int c) { return c >= 'a' && c <= 'z'; }
    bool isUppercase(int c) { return c >= 'A' && c <= 'Z'; }
    bool isDigit(int c) { return c >= '0' && c <= '9'; }
    bool isLetterOrDigit(int c) { return isLowercase(c) || isUppercase(c)  || isDigit(c); }
}

////////////////////////////////////////////////////////////////////

void NameManager::insert(string name)
{
    for (auto c : name) if (!isLetterOrDigit(c)) throw FATALERROR("Name can contain only letters and digits");

    if (isUppercase(name[0])) _globalSet.insert(name);
    else if (isLowercase(name[0])) _localSetStack.top().insert(name);
    else throw FATALERROR("First character in name must be a letter");
}

////////////////////////////////////////////////////////////////////

bool NameManager::evaluateBoolean(string expression)
{
    return BooleanExpression::evaluateBoolean(expression, [this](string name)
    {
        return _globalSet.count(name) > 0 || _localSetStack.top().count(name) > 0;
    });
}

////////////////////////////////////////////////////////////////////

string NameManager::evaluateConditionalValue(string expression)
{
    return BooleanExpression::evaluateConditionalValue(expression, [this](string name)
    {
        return _globalSet.count(name) > 0 || _localSetStack.top().count(name) > 0;
    });
}

////////////////////////////////////////////////////////////////////
