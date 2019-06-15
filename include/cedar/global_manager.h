#ifndef CEDAR_GLOBAL_MANAGER_H
#define CEDAR_GLOBAL_MANAGER_H

#include <tuple>
#include <memory>

#include <cedar/cxx14.h>
#include <cedar/type_list.h>
#include <cedar/global_params.h>
#include <cedar/util/basic_logger.h>
#include <cedar/util/basic_allocator.h>
#include <cedar/util/agg_timer.h>

namespace cedar {

template<class U> class global_manager;
template<typename... types>
class global_manager<type_list<types...>>
{
public:
	using gtype = std::tuple<std::unique_ptr<types>...>;
	using logger = typename std::tuple_element<0, gtype>::type::element_type;
	using timer  = typename std::tuple_element<1, gtype>::type::element_type;
	using memory = typename std::tuple_element<2, gtype>::type::element_type;
	global_manager() {}
	global_manager(global_params & params) :
		globals(std::make_unique<types>(params)...) {}
	template<class T>
	T & get()
	{
		return *(std::get<std::unique_ptr<T>>(globals));
	}

protected:
	gtype globals;
};

using reg_globals = type_list<basic_logger, agg_timer, basic_allocator>;
using gmant = global_manager<reg_globals>;
extern global_manager<reg_globals> gman;

}

#endif
