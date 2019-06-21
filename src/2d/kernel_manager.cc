#include <cedar/2d/relax.h>
#include <cedar/2d/coarsen.h>
#include <cedar/2d/interp.h>
#include <cedar/2d/residual.h>
#include <cedar/2d/restrict.h>
#include <cedar/2d/solve_cg.h>
#include <cedar/2d/kernel_manager.h>
#include <cedar/malloc_pool.h>

namespace cedar { namespace cdr2 {
kman_ptr build_kernel_manager(config & conf)
{
	return build_kernel_manager(build_kernel_params(conf));
}


kman_ptr build_kernel_manager(std::shared_ptr<kernel_params> params)
{
	log::status << *params << std::endl;

	auto kman = std::make_shared<kernel_manager<kernels2,stypes>>(params);
	kman->add<point_relax, rbgs>("system", kmode::serial);
	kman->add<point_relax, rbgs>("offload", kmode::offload);
	kman->add<point_relax, rbgs>("omp", kmode::omp);
	kman->add<line_relax<relax_dir::x>, lines<relax_dir::x>>("system");
	kman->add<line_relax<relax_dir::y>, lines<relax_dir::y>>("system");
	kman->add<coarsen_op, galerkin>("system");
	kman->add<interp_add, interp_f90>("system", kmode::serial);
	kman->add<interp_add, interp_f90>("offload", kmode::offload);
	kman->add<interp_add, interp_f90>("omp", kmode::omp);
	kman->add<restriction, restrict_f90>("system", kmode::serial);
	kman->add<restriction, restrict_f90>("offload", kmode::offload);
	kman->add<restriction, restrict_f90>("omp", kmode::omp);
	kman->add<residual, residual_f90>("system", kmode::serial);
	kman->add<residual, residual_f90>("offload", kmode::offload);
	kman->add<residual, residual_f90>("omp", kmode::omp);
	kman->add<setup_interp, setup_interp_f90>("system");
	kman->add<solve_cg, solve_cg_f90>("system", kmode::serial);
	kman->add<solve_cg, solve_cg_f90>("offload", kmode::offload);

	kman->set<point_relax>("system");
	kman->set<line_relax<relax_dir::x>>("system");
	kman->set<line_relax<relax_dir::y>>("system");
	kman->set<coarsen_op>("system");
	kman->set<interp_add>("system");
	kman->set<restriction>("system");
	kman->set<residual>("system");
	kman->set<setup_interp>("system");
	kman->set<solve_cg>("system");

	if (params->offload) {
		kman->set<residual>("offload");
		kman->set<restriction>("offload");
		kman->set<point_relax>("offload");
		kman->set<interp_add>("offload");
		kman->set<solve_cg>("offload");
	}

	// register services
	auto & services = kman->services();
	services.add<mempool, malloc_pool>("malloc");
	services.set<mempool>("malloc");

	return kman;
}
}}
