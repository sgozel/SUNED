// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SG_NUMA_H
#define SG_NUMA_H

#include <vector>

#ifdef SG_USE_NUMA

#include <new>

// Minimal uninitialized allocator — avoids value-init touching pages
template<typename T>
struct uninit_allocator : std::allocator<T> {
    template<typename U> struct rebind { using other = uninit_allocator<U>; };
    void construct(T*) {} // no-op: skip zero-initialization
    template<typename... Args>
    void construct(T* p, Args&&... args) {
        ::new((void*)p) T(std::forward<Args>(args)...);
    }
};


template<typename T>
using sg_vec = std::vector<T, uninit_allocator<T>>;

#else

template<typename T>
using sg_vec = std::vector<T>;

#endif

#endif
