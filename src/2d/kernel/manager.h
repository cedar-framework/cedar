#ifndef BOXMG_2D_KERNEL_MANAGER_H
#define BOXMG_2D_KERNEL_MANAGER_H

#include <memory>

#include "boxmg-common.h"


namespace boxmg { namespace bmg2d { namespace kernel {

class Manager
{
public:
	template <typename S, typename T>
		static void add(S&& kname, const T&cmd);
	template <typename S, typename T>
		static void reg(const std::string & kname, S&& kid, const T&cmd);
	template <typename T>
		static void set(const std::string & kname, const T& cmd);
	static void change(const std::string & kname, const std::string & rname);
	template <typename... Args>
		static void run(const std::string & kname, Args&&... args);
private:
		static void check_init();
		static void init_reg();
		static void init();
		static std::unique_ptr<boxmg::KernelManager> instance;
		static std::unique_ptr<std::map<std::string,boxmg::KernelManager>> avail;

};


template<typename S, typename T>
void Manager::add(S&& kname, const T&cmd)
{
	check_init();
	instance->add(std::forward<S>(kname), cmd);
}


template<typename T>
void Manager::set(const std::string & kname, const T & cmd)
{
	check_init();
	instance->set(kname, cmd);
}


template<typename S, typename T>
	void Manager::reg(const std::string & kname, S&& kid, const T&cmd)
{
	if (!avail)
		init_reg();

	(*avail)[kname].add(std::forward<S>(kid), cmd);
}


template <typename... Args>
	void Manager::run(const std::string & kname, Args&&... args)
{
	check_init();
	instance->run(kname, std::forward<decltype(args)>(args)...);
}

}}}

#endif
