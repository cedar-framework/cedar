#ifndef CEDAR_LINE_RELAX_H
#define CEDAR_LINE_RELAX_H

#include <cedar/types.h>
#include <cedar/kernel.h>

namespace cedar { namespace kernels {

	template<class solver_types, relax_dir rdir>
	class line_relax : public kernel<solver_types>
	{
	public:
		template<class sten>
		using stencil_op = typename kernel<solver_types>::template stencil_op<sten>;
		using comp_sten = typename kernel<solver_types>::comp_sten;
		using full_sten = typename kernel<solver_types>::full_sten;
		using grid_func = typename kernel<solver_types>::grid_func;
		using relax_stencil = typename kernel<solver_types>::relax_stencil;

		const static std::string name() { return "line relaxation"; }

		virtual void setup(const stencil_op<comp_sten> & so,
		                   relax_stencil & sor) = 0;
		virtual void setup(const stencil_op<full_sten> & so,
		                   relax_stencil & sor) = 0;


		virtual void run(const stencil_op<comp_sten> & so,
		                 grid_func & x,
		                 const grid_func & b,
		                 const relax_stencil & sor,
		                 grid_func & res,
		                 cycle::Dir cdir) = 0;
		virtual void run(const stencil_op<full_sten> & so,
		                 grid_func & x,
		                 const grid_func & b,
		                 const relax_stencil & sor,
		                 grid_func & res,
		                 cycle::Dir cdir) = 0;
	};
}}

#endif
