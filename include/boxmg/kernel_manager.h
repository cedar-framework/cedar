#ifndef BOXMG_KERNEL_MANAGER_H
#define BOXMG_KERNEL_MANAGER_H

#include <memory>
#include <string>
#include <map>
#include <exception>

#include "kernel.h"
#include "util/timer.h"
#include "util/log.h"

namespace boxmg {

class kernel_manager
{
public:
	template <typename S, typename T>
		void add(S&& kname, const T& cmd)
	{
		kmap.emplace(std::make_pair(std::forward<S>(kname), std::make_shared<T>(cmd)));
	}

	template <typename T>
		void set(const std::string & kname, const T & cmd)
	{
		kmap_t::const_iterator it = kmap.find(kname);
		if (it != kmap.end()) {
			T *c = dynamic_cast<T*>(it->second.get());
			if (c) {
				kmap[kname] = std::make_shared<T>(cmd);
			} else {
				throw std::invalid_argument("Argument does not match signature of kernel: " + kname);
			}
		} else {
			throw std::invalid_argument("kernel not found: " + kname);
		}
	}

	template <typename... Args>
		void run(const std::string & kname, Args&&... args)
	{
		kmap_t::const_iterator it = kmap.find(kname);
		if (it != kmap.end()) {
			kernel<Args...> *c = dynamic_cast<kernel<Args...>*>(it->second.get());
			if (c) {
				log::debug << "Running kernel: " << kname << std::endl;
				Timer kern_timer(kname);
				if (log::debug.active())
					kern_timer.begin();
				(*c)(std::forward<decltype(args)>(args)...);
				if (log::debug.active())
					kern_timer.end();
			} else {
				std::string msg("Incorrect arguments for kernel: " + kname);
				log::error << msg << std::endl;
				//throw std::invalid_argument(msg);
			}
		} else {
			log::error << "kernel not found: "  << kname << std::endl;
			//throw std::invalid_argument("kernel not found: " + kname);
		}
	}

	std::shared_ptr<kernel_base> & operator[](const std::string & kname)
	{
		return kmap[kname];
	}


	std::shared_ptr<kernel_base> at(const std::string & kname)
	{
		try {
			return kmap.at(kname);
		} catch(const std::out_of_range  &ex) {
			log::error << "Cound not find kernel: " << kname << std::endl;
			return kmap[kname];
		}
	}
private:
	using kmap_t = std::map<std::string, std::shared_ptr<kernel_base>>;
	kmap_t kmap;
};

}

#endif
