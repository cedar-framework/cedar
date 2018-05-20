#ifndef CEDAR_INTERP_ADD_H
#define CEDAR_INTERP_ADD_H

#include <cedar/kernel.h>

namespace cedar { namespace kernels {

	template<class solver_types>
	class interp_add : public kernel<solver_types>
	{
	public:
		using grid_func = typename kernel<solver_types>::grid_func;
		using prolong_op = typename kernel<solver_types>::prolong_op;

		const static std::string name() { return "interpolation and add"; }

		virtual void run(const prolong_op & P,
		                 const grid_func & coarse,
		                 const grid_func & residual,
		                 grid_func & fine) = 0;
	};
}}

#endif
