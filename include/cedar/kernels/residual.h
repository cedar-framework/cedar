#ifndef CEDAR_RESIDUAL_H
#define CEDAR_RESIDUAL_H

#include <cedar/kernel.h>

namespace cedar { namespace kernels {
	template<class solver_types>
	class residual : public kernel<solver_types>
	{
	public:
		template<class sten>
		using stencil_op = typename kernel<solver_types>::template stencil_op<sten>;
		using comp_sten = typename kernel<solver_types>::comp_sten;
		using full_sten = typename kernel<solver_types>::full_sten;
		using grid_func = typename kernel<solver_types>::grid_func;

		const std::string name = "residual";

		virtual void run(const stencil_op<comp_sten> & so,
		                 const grid_func & x,
		                 const grid_func & b,
		                 grid_func & r) = 0;
		virtual void run(const stencil_op<full_sten> & so,
		                 const grid_func & x,
		                 const grid_func & b,
		                 grid_func & r) = 0;
	};
}}

#endif
