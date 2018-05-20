#ifndef CEDAR_POINT_RELAX_H
#define CEDAR_POINT_RELAX_H

#include <cedar/kernel.h>
#include <cedar/types.h>

namespace cedar { namespace kernels {

	template<class solver_types>
	class point_relax : public kernel<solver_types>
	{
	public:
		template<class sten>
		using stencil_op = typename kernel<solver_types>::template stencil_op<sten>;
		using comp_sten = typename kernel<solver_types>::comp_sten;
		using full_sten = typename kernel<solver_types>::full_sten;
		using relax_stencil = typename kernel<solver_types>::relax_stencil;
		using grid_func = typename kernel<solver_types>::grid_func;

		const static std::string name() { return "point relaxation"; }

		point_relax() {}
		point_relax(std::function<void(const stencil_op<comp_sten>&, relax_stencil &)> scomp,
		            std::function<void(const stencil_op<full_sten>&, relax_stencil &)> sfull,
		            std::function<void(const stencil_op<comp_sten>&, grid_func &, const grid_func &,
		                               const relax_stencil &, cycle::Dir)> rcomp,
		            std::function<void(const stencil_op<comp_sten>&, grid_func &, const grid_func &,
		                               const relax_stencil &, cycle::Dir)> rfull) :
			setup_comp(scomp), setup_full(sfull), run_comp(rcomp), run_full(rfull) {}

		virtual void setup(const stencil_op<comp_sten> & so,
		                   relax_stencil & sor)
		{
			if (setup_comp)
				setup_comp(so, sor);
			else
				log::error << name() << ": setup routine not provided" << std::endl;

		}
		virtual void setup(const stencil_op<full_sten> & so,
		                   relax_stencil & sor)
		{
			if (setup_full)
				setup_full(so, sor);
			else
				log::error << name() << ": setup routine not provided" << std::endl;
		}


		virtual void run(const stencil_op<comp_sten> & so,
		                 grid_func & x,
		                 const grid_func & b,
		                 const relax_stencil & sor,
		                 cycle::Dir cdir)
		{
			if (run_comp)
				run_comp(so, x, b, sor, cdir);
			else
				log::error << name() << ": run routine not provided" << std::endl;
		}
		virtual void run(const stencil_op<full_sten> & so,
		                 grid_func & x,
		                 const grid_func & b,
		                 const relax_stencil & sor,
		                 cycle::Dir cdir)
		{
			if (run_full)
				run_full(so, x, b, sor, cdir);
			else
				log::error << name() << ": run routine not provided" << std::endl;
		}

	protected:
		std::function<void(const stencil_op<comp_sten>&, relax_stencil &)> setup_comp;
		std::function<void(const stencil_op<full_sten>&, relax_stencil &)> setup_full;
		std::function<void(const stencil_op<comp_sten>&, grid_func &, const grid_func &,
		                   const relax_stencil &, cycle::Dir)> run_comp;
		std::function<void(const stencil_op<full_sten>&, grid_func &, const grid_func &,
		                   const relax_stencil &, cycle::Dir)> run_full;
	};
}}
#endif
