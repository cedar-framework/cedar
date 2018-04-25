#ifndef CEDAR_3D_REDIST_CHOLESKY_SOLVER_H
#define CEDAR_3D_REDIST_CHOLESKY_SOLVER_H

#include <memory>

#include <cedar/config.h>
#include <cedar/3d/types.h>
#include <cedar/3d/kernel_manager.h>

namespace cedar { namespace cdr3 {

class cholesky_solver
{
public:
	static const bool is_serial = true;
	using stypes = cdr3::stypes;
	using stencil_op = typename stypes::template stencil_op<xxvii_pt>;
	using grid_func = typename stypes::grid_func;
	cholesky_solver(stencil_op & sop, std::shared_ptr<config> conf);
	void cycle(grid_func & x, const grid_func & b);
	config & get_config() { return *conf; }

protected:
	kman_ptr kman;
	grid_func ABD;
	real_t *bbd;
	std::shared_ptr<config> conf;
};

}}

#endif
