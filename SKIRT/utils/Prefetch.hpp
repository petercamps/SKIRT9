/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PREFETCH_HPP
#define PREFETCH_HPP

#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** This namespace offers functions that ask the processor to load memory into its cache in the
    background, so that the memory is (more likely to be) available by the time it is actually
    accessed. This can hide the latency of memory accesses that cannot be predicted by the
    processor itself, such as following links between the nodes of a large tree that does not fit
    in the cache.

    A prefetch request is just a hint without any other effect: it does not wait for the memory to
    arrive, and it can never fail, even for an invalid address. Prefetching is not part of standard
    C++. The builtin function offering it is available in GCC and Clang, and in compilers that
    identify as one of these (such as the Intel compilers on Linux and macOS, and the Clang-based
    Intel oneAPI compilers on all platforms); for other compilers, including MSVC, the functions in
    this namespace do nothing.

    Note: conditional compilation for system-dependent code should normally be restricted to the
    System class. This header is a deliberate exception, because the functions must be inlined at
    the call site to be useful. */
namespace Prefetch
{
    /** This function asks the processor to load the cache line containing the specified address
        into its cache in the background. */
#if defined(__GNUC__) || defined(__clang__)
    inline void address(const void* address)
    {
        __builtin_prefetch(address);
    }
#else
    inline void address(const void* /*address*/) {}
#endif

    /** This function asks the processor to load all cache lines holding the specified object into
        its cache in the background. The function assumes cache lines of at least 64 bytes; it
        requests the first byte of the object, every 64th byte after it, and the last byte, so that
        all cache lines are covered regardless of the alignment of the object. For a small object,
        the compiler can unroll the loop completely, because the size of the object is known at
        compile time. */
    template<class T> inline void object(const T* object)
    {
        const char* first = static_cast<const char*>(static_cast<const void*>(object));
        for (size_t offset = 0; offset < sizeof(T); offset += 64) address(first + offset);
        address(first + sizeof(T) - 1);
    }
}

////////////////////////////////////////////////////////////////////

#endif
