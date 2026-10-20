/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPAN_HPP
#define SPAN_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** A Span instance refers to a contiguous sequence of values of type \c T that it does not own. It
    offers a minimal subset of the \c std::span class template, which requires C++20. Copying a
    Span is cheap: it copies the reference to the values, not the values themselves. The values
    can be modified through the Span unless \c T is a const type. A Span remains valid only as long
    as the values it refers to are not moved or released by their owner.

    A Span on non-constant values converts implicitly to a Span on the same values declared
    constant, so that a writable view can be passed where a read-only view is expected. The begin()
    and end() functions allow using a Span with the standard algorithms and in a range-based for
    loop. */
template<typename T> class Span
{
public:
    /** This constructor creates a Span referring to the specified number of values starting at the
        specified address. */
    Span(T* first, size_t size) : _first{first}, _size{size} {}

    /** This constructor converts a Span on non-constant values into a Span on the same values
        declared constant. It is available only if \c T is a const type. */
    template<typename U, typename = std::enable_if_t<std::is_same<const U, T>::value && !std::is_same<U, T>::value>>
    Span(const Span<U>& other) : _first{other.begin()}, _size{other.size()}
    {}

    /** This function returns the number of values in the Span. */
    size_t size() const { return _size; }

    /** This function returns a reference to the value with the specified index, which must be
        smaller than the number of values in the Span. */
    T& operator[](size_t i) const { return _first[i]; }

    /** This function returns a pointer to the first value in the Span. */
    T* begin() const { return _first; }

    /** This function returns a pointer just beyond the last value in the Span. */
    T* end() const { return _first + _size; }

private:
    T* _first{nullptr};
    size_t _size{0};
};

////////////////////////////////////////////////////////////////////

#endif
