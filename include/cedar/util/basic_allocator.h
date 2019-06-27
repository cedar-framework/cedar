#ifndef CEDAR_UTIL_BASIC_ALLOCATOR_H
#define CEDAR_UTIL_BASIC_ALLOCATOR_H

#include <unordered_map>
#include <functional>
#include <cedar/global_params.h>
#include <cedar/memory_types.h>

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
	template<class T>
	memory::location hint(const T *addr)
	{
		return hint_impl(static_cast<const void*>(addr));
	}
	std::size_t prefetchsize() { return pfsize; }
	template<class T>
	void save(const std::string & key, T * value)
	{
		saved[key] = static_cast<void*>(value);
	}
	template<class T>
	T * load(const std::string & key)
	{
		auto it = saved.find(key);
		if (it != saved.end())
			return static_cast<T*>(it->second);
		else
			return static_cast<T*>(nullptr);
	}
protected:
	std::size_t pfsize;
	std::function<void*(std::size_t nbytes)> allocator;
	std::function<void(void*addr)> deallocator;
	std::unordered_map<const void*, memory::location> hints;
	std::unordered_map<std::string, void*> saved;

	void prefetch_impl(void *addr, std::size_t nbytes);
	memory::location hint_impl(const void *addr);
};

}

#endif
