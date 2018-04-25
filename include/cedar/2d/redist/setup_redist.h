#ifndef CEDAR_2D_REDIST_SETUP_REDIST_H
#define CEDAR_2D_REDIST_SETUP_REDIST_H

#include <cedar/kernel_params.h>
#include <cedar/2d/mpi/kernel_manager.h>
#include <cedar/2d/mpi/redist_solver.h>

namespace cedar { namespace cdr2 { namespace mpi {

template<class inner_solver>
std::function<void(mpi::grid_func &, const mpi::grid_func &)>
	create_redist_solver(kman_ptr kman,
	                     config & conf,
	                     mpi::stencil_op<nine_pt> & cop,
	                     std::shared_ptr<config> cg_conf,
	                     std::array<int, 2> & choice)
{
	auto params = build_kernel_params(conf);
	using rsolver = mpi::redist_solver<inner_solver>;
	auto cg_bmg = std::make_shared<rsolver>(cop,
	                                        kman->get_ptr<halo_exchange>().get(),
	                                        cg_conf,
	                                        choice);

	auto coarse_solver = [=](mpi::grid_func & x, const mpi::grid_func & b)
	{
		cg_bmg->solve(x, b);
		if (params->per_mask())
			kman->run<halo_exchange>(x);
	};

	return coarse_solver;
}

}}}

#endif
