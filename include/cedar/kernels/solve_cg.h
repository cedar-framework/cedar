#ifndef CEDAR_SOLVE_CG
#define CEDAR_SOLVE_CG

#include <cedar/kernel.h>

namespace cedar { namespace kernels {

	template<class solver_types>
	class solve_cg : public kernel<solver_types>
	{
	public:
		template<class sten>
		using stencil_op = typename kernel<solver_types>::template stencil_op<sten>;
		using comp_sten = typename kernel<solver_types>::comp_sten;
		using full_sten = typename kernel<solver_types>::full_sten;
		using grid_func = typename kernel<solver_types>::grid_func;

		const static std::string name() { return "coarse-grid solve"; }

		virtual void setup(const stencil_op<comp_sten> & so,
		                   grid_func & ABD) = 0;
		virtual void setup(const stencil_op<full_sten> & so,
		                   grid_func & ABD) = 0;

		virtual void run(grid_func & x,
		                 const grid_func & b,
		                 const grid_func &ABD,
		                 real_t *bbd) = 0;
	};
}}

#endif
