#ifndef CEDAR_KERNEL_MANAGER_H
#define CEDAR_KERNEL_MANAGER_H

#include <cedar/type_list.h>
#include <cedar/kernel_params.h>
#include <cedar/config.h>


namespace cedar {

template<class tlist>
class kernel_manager : public type_map<tlist>
{
public:
	using parent = type_map<tlist>;

    kernel_manager(std::shared_ptr<kernel_params> params) : params(params) {}
	kernel_manager(config & conf) { params = build_kernel_params(conf); }

	template<class T, class rclass>
	void add(const std::string & name)
	{
		parent::template add<T, rclass>(name);
		this->init<T>(name);
	}


	template<class T, class rclass, class... Args>
	void add(const std::string & name, Args&&... args)
	{
		parent::template add<T, rclass>(name, std::forward<Args>(args)...);
		this->init<T>(name);
	}


	template<class T, class... Args>
	void setup(Args&&... args)
	{
		auto & kern = parent::template get<T>();
		log::debug << "setup kernel <" << kern.name << ">" << std::endl;
		kern.setup(std::forward<Args>(args)...);
	}

	template<class T, class... Args>
	void run(Args&&... args)
	{
		auto & kern = parent::template get<T>();
		log::debug << "running kernel <" << kern.name << ">" << std::endl;
		kern.run(std::forward<Args>(args)...);
	}


	template<class T>
	void init(const std::string & name)
	{
		parent::template get<T>(name).add_params(params);
	}

protected:
	std::shared_ptr<kernel_params> params;
};

}

#endif
