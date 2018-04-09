#include <cedar/3d/relax.h>
#include <cedar/3d/coarsen.h>
#include <cedar/3d/interp.h>
#include <cedar/3d/residual.h>
#include <cedar/3d/restrict.h>
#include <cedar/3d/solve_cg.h>
#include <cedar/3d/kernel_manager.h>

namespace cedar { namespace cdr3 {
kman_ptr build_kernel_manager(config::reader & conf)
{
	return build_kernel_manager(build_kernel_params(conf));
}


kman_ptr build_kernel_manager(std::shared_ptr<kernel_params> params)
{
	auto kman = std::make_unique<kernel_manager<kernels3>>(params);
	kman->add<point_relax, rbgs>("system");
	kman->add<coarsen_op, galerkin>("system");
	kman->add<interp_add, interp_f90>("system");
	kman->add<restriction, restrict_f90>("system");
	kman->add<residual, residual_f90>("system");
	kman->add<setup_interp, setup_interp_f90>("system");
	kman->add<solve_cg, solve_cg_f90>("system");

	kman->set<point_relax>("system");
	kman->set<coarsen_op>("system");
	kman->set<interp_add>("system");
	kman->set<restriction>("system");
	kman->set<residual>("system");
	kman->set<setup_interp>("system");
	kman->set<solve_cg>("system");

	return kman;
}
}}
