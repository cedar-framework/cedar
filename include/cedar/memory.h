#ifndef CEDAR_MEMORY_H
#define CEDAR_MEMORY_H

#include <functional>
#include <cedar/config.h>
namespace cedar { namespace memory {

// TODO: move these to global manager
extern std::function<void*(std::size_t nbytes)> allocator;
extern std::function<void(void * addr)> deallocator;

void init(config & conf);
template<class T>
inline T * alloc(std::size_t N) { return static_cast<T*>(allocator(N*sizeof(T))); }
template<class T>
inline void free(T *addr) { deallocator(static_cast<void*>(addr)); }

}}

#endif
