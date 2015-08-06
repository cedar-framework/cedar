#ifndef BOXMG_KERNEL_REGISTRY_H
#define BOXMG_KERNEL_REGISTRY_H

#include <memory>
#include <map>
#include "kernel_manager.h"

namespace boxmg {


class KernelRegistry
{
public:
	template <typename S, typename T>
		void add(const std::string & kname, S&& kid, const T&cmd);
	void set(const std::string & kname, const std::string & kid);
	template <typename... Args>
		void run(const std::string & kname, Args&&... args);
protected:
	KernelManager active;
	std::map<std::string, KernelManager> avail;
};


template<typename S, typename T>
void KernelRegistry::add(const std::string & kname, S&& kid, const T&cmd)
{
	avail[kname].add(std::forward<S>(kid), cmd);
}


template <typename... Args>
	void KernelRegistry::run(const std::string & kname, Args&&... args)
{
	active.run(kname, std::forward<decltype(args)>(args)...);
}

}

#endif
