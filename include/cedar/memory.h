#ifndef CEDAR_MEMORY_H
#define CEDAR_MEMORY_H

#include <cedar/global_manager.h>
#include <cedar/memory_types.h>

namespace cedar { namespace memory {

using mtype = global_manager<reg_globals>::memory;

template<class T>
inline T * alloc(std::size_t N) { return gman.get<mtype>().alloc<T>(N); }

template<class T>
inline void free(T *addr) { return gman.get<mtype>().free(addr); }

template<class T>
inline void prefetch(T *addr, std::size_t N) { gman.get<mtype>().prefetch(addr, N); }

inline void sync() { gman.get<mtype>().sync(); }

template<class T>
inline location hint(const T *addr) { return gman.get<mtype>().hint(addr); }

inline std::size_t prefetchsize() { return gman.get<mtype>().prefetchsize(); }
}}

#endif
