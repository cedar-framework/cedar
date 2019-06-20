#ifndef CEDAR_KERNEL_H
#define CEDAR_KERNEL_H

#include <memory>
#include <cedar/kernel_params.h>
#include <cedar/service_manager.h>

namespace cedar {

enum class kmode { serial, omp, offload };

	template<class solver_types>
	class kernel
	{
	public:
		// extract types from solver_types
		template<class sten>
			using stencil_op = typename solver_types::template stencil_op<sten>;
		using comp_sten = typename solver_types::comp_sten;
		using full_sten = typename solver_types::full_sten;
		using grid_func = typename solver_types::grid_func;
		using prolong_op = typename solver_types::prolong_op;
		using restrict_op = typename solver_types::restrict_op;
		using relax_stencil = typename solver_types::relax_stencil;


		void add_params(std::shared_ptr<kernel_params> params)
		{
			this->params = params;
		}

		void add_services(service_manager<solver_types> *services)
		{
			this->services = services;
		}

	protected:
		std::shared_ptr<kernel_params> params;
		service_manager<solver_types> *services;
	};
}

#endif
