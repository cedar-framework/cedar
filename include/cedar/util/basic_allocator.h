#ifndef CEDAR_UTIL_BASIC_ALLOCATOR_H
#define CEDAR_UTIL_BASIC_ALLOCATOR_H

#include <functional>
#include <cedar/global_params.h>

namespace cedar
{

class basic_allocator
{
public:
	basic_allocator(global_params & params);
	template<class T>
	T * alloc(std::size_t N) { return static_cast<T*>(allocator(N*sizeof(T))); }
	template<class T>
	void free(T *addr) { deallocator(static_cast<void*>(addr)); }
	template<class T>
	void prefetch(T *addr, std::size_t n)
	{
		prefetch_impl(static_cast<void*>(addr), n*sizeof(T));
	}
	void sync();
protected:
	std::function<void*(std::size_t nbytes)> allocator;
	std::function<void(void*addr)> deallocator;
	void prefetch_impl(void *addr, std::size_t nbytes);
};

}

#endif
