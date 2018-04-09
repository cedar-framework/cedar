#ifndef CEDAR_MATVEC_H
#define CEDAR_MATVEC_H

#include <cedar/kernel.h>

namespace cedar { namespace kernels {

	template<class solver_types>
	class matvec : public kernel<solver_types>
	{
	public:
		template<class sten>
		using stencil_op = typename kernel<solver_types>::template stencil_op<sten>;
		using comp_sten = typename kernel<solver_types>::comp_sten;
		using full_sten = typename kernel<solver_types>::full_sten;
		using grid_func = typename kernel<solver_types>::grid_func;

		const std::string name = "matvec";

		virtual void run(const stencil_op<comp_sten> & so,
		                 const grid_func & x,
		                 grid_func & y)
		{
			if (run_comp)
				run_comp(so, x, y);
			else
				log::error << name << ": routine not provided" << std::endl;
		}
		virtual void run(const stencil_op<full_sten> & so,
		                 const grid_func & x,
		                 grid_func & y)
		{
			if (run_full)
				run_full(so, x, y);
			else
				log::error << name << ": routine not provided" << std::endl;
		}

	protected:
		std::function<void(const stencil_op<comp_sten> &, const grid_func &, grid_func &)> run_comp;
		std::function<void(const stencil_op<full_sten> &, const grid_func &, grid_func &)> run_full;
	};
}}

#endif

