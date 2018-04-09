#ifndef CEDAR_SETUP_INTERP_H
#define CEDAR_SETUP_INTERP_H

#include <cedar/kernel.h>


namespace cedar { namespace kernels {

	template<class solver_types>
	class setup_interp : public kernel<solver_types>
	{
	public:
		template<class sten>
		using stencil_op = typename kernel<solver_types>::template stencil_op<sten>;
		using comp_sten = typename kernel<solver_types>::comp_sten;
		using full_sten = typename kernel<solver_types>::full_sten;
		using prolong_op = typename kernel<solver_types>::prolong_op;

		const std::string name = "setup interpolation";

		setup_interp() {}
		setup_interp(std::function<void(const stencil_op<comp_sten> &,
		                                const stencil_op<full_sten> &,
		                                prolong_op&)> crun,
		             std::function<void(const stencil_op<full_sten> &,
		                                const stencil_op<full_sten> &,
		                                prolong_op&)> frun) :
			run_comp(crun), run_full(frun) {}

		virtual void run(const stencil_op<comp_sten> & fop,
		                 const stencil_op<full_sten> & cop,
		                 prolong_op & P)
		{
			if (run_comp)
				run_comp(fop, cop, P);
			else
				log::error << name << ": routine not provided" << std::endl;
		}
		virtual void run(const stencil_op<full_sten> & fop,
		                 const stencil_op<full_sten> & cop,
		                 prolong_op & P)
		{
			if (run_full)
				run_full(fop, cop, P);
			else
				log::error << name << ": routine not provided" << std::endl;
		}

	protected:
		std::function<void(const stencil_op<comp_sten> &,
		                   const stencil_op<full_sten> &,
		                   prolong_op&)> run_comp;
		std::function<void(const stencil_op<full_sten> &,
		                   const stencil_op<full_sten> &,
		                   prolong_op&)> run_full;
	};
}}

#endif
