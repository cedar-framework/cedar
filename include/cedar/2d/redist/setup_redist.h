#ifndef CEDAR_2D_REDIST_SETUP_REDIST_H
#define CEDAR_2D_REDIST_SETUP_REDIST_H

#include <cedar/kernel_params.h>
#include <cedar/2d/mpi/kernel_manager.h>
#include <cedar/2d/mpi/redist_solver.h>

namespace cedar { namespace cdr2 { namespace mpi {

template<class inner_solver>
std::shared_ptr<mpi::redist_solver<inner_solver>> create_redist_ptr(kman_ptr kman,
                                                                    mpi::stencil_op<nine_pt> & cop,
                                                                    std::shared_ptr<config> cg_conf,
                                                                    std::array<int, 2> & choice)
{

	using rsolver = mpi::redist_solver<inner_solver>;
	auto cg_bmg = std::make_shared<rsolver>(cop,
	                                        kman->services_ptr(),
	                                        cg_conf,
	                                        choice);
	return cg_bmg;
}


template<class inner_solver>
std::function<void(mpi::grid_func&, const mpi::grid_func &)>
wrap_redist_ptr(config & conf, kman_ptr kman, std::shared_ptr<mpi::redist_solver<inner_solver>> cg_bmg)
{
	auto params = build_kernel_params(conf);

	auto coarse_solver = [=](mpi::grid_func & x, const mpi::grid_func & b)
	{
		cg_bmg->solve(x, b);
		if (params->per_mask()) {
			auto & halo_service = kman->services().get<halo_exchange>();
			halo_service.run(x);
		}
	};

	return coarse_solver;
}


template<class inner_solver>
std::function<void(mpi::grid_func &, const mpi::grid_func &)>
	create_redist_solver(kman_ptr kman,
	                     config & conf,
	                     mpi::stencil_op<nine_pt> & cop,
	                     std::shared_ptr<config> cg_conf,
	                     std::array<int, 2> & choice)
{
	auto cg_bmg = create_redist_ptr<inner_solver>(kman, cop, cg_conf, choice);

	return wrap_redist_ptr<inner_solver>(conf, kman, cg_bmg);
}


}}}

#endif
