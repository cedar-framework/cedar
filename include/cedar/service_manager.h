#ifndef CEDAR_SERVICE_MANAGER_H
#define CEDAR_SERVICE_MANAGER_H

#include <cedar/type_list.h>
#include <cedar/log.h>
#include <cedar/services/halo_exchange.h>
#include <cedar/services/message_passing.h>
#include <cedar/services/mempool.h>

namespace cedar {

template<class tlist>
class service_manager_gen : public type_map<tlist>
{
public:
	using parent = type_map<tlist>;
	using mtype = typename parent::mtype;
	using parent::kerns;

	service_manager_gen(std::shared_ptr<kernel_params> params) : params(params) {}

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

	template<class T>
	void init(const std::string & name)
	{
		auto & kern = parent::template get<T>(name);
		kern.add_params(params);
	}

	template<class T>
	T* fortran_handle()
	{
		auto kern = parent::template get_ptr<T>();
		if (!kern)
			log::error << "service not found: " << T::name() << std::endl;
		return kern.get();
	}


	void set_user_reg(std::function<void(service_manager_gen<tlist> &)> reg) {
		if (reg)
			reg(*this);
		user_reg = reg;
	}


	std::function<void(service_manager_gen<tlist> &)> get_user_reg() {
		return user_reg;
	}


protected:
	std::shared_ptr<kernel_params> params;
	std::function<void(service_manager_gen<tlist> &)> user_reg;
};


template<class solver_types>
using reg_services = type_list<services::halo_exchange<solver_types>,
                               services::message_passing,
                               services::mempool>;

template<class solver_types>
using service_manager = service_manager_gen<reg_services<solver_types>>;

}
#endif
