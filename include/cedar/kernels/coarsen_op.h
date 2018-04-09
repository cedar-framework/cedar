#ifndef CEDAR_COARSEN_OP_H
#define CEDAR_COARSEN_OP_H

namespace cedar { namespace kernels {
	template<class solver_types>
	class coarsen_op : public kernel<solver_types>
	{
	public:
		template<class sten>
		using stencil_op = typename kernel<solver_types>::template stencil_op<sten>;
		using comp_sten = typename kernel<solver_types>::comp_sten;
		using full_sten = typename kernel<solver_types>::full_sten;
		using grid_func = typename kernel<solver_types>::grid_func;
		using prolong_op = typename kernel<solver_types>::prolong_op;

		coarsen_op() {}
		coarsen_op(std::function<void(const prolong_op &, const stencil_op<comp_sten>, stencil_op<full_sten>)> rcomp,
		           std::function<void(const prolong_op &, const stencil_op<full_sten>, stencil_op<full_sten>)> rfull):
			run_comp(rcomp), run_full(rfull) {}

		const std::string name = "coarsen operator";

		virtual void run(const prolong_op & P,
		                  const stencil_op<comp_sten> & fop,
		                  stencil_op<full_sten> & cop)
		{
			if (run_comp)
				run_comp(P, fop, cop);
			else
				log::error << name << ": routine not provided" << std::endl;
		}
		virtual void run(const prolong_op & P,
		                  const stencil_op<full_sten> & fop,
		                  stencil_op<full_sten> & cop)
		{
			if (run_full)
				run_full(P, fop, cop);
			else
				log::error << name << ": routine not provided" << std::endl;
		}

	protected:
		std::function<void(const prolong_op &, const stencil_op<comp_sten>, stencil_op<full_sten>)> run_comp;
		std::function<void(const prolong_op &, const stencil_op<full_sten>, stencil_op<full_sten>)> run_full;
	};
}}

#endif
