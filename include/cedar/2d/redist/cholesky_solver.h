#ifndef CEDAR_2D_REDIST_CHOLESKY_SOLVER_H
#define CEDAR_2D_REDIST_CHOLESKY_SOLVER_H

#include <memory>

#include <cedar/2d/types.h>
#include <cedar/config.h>
#include <cedar/2d/kernel_manager.h>

namespace cedar { namespace cdr2 {

class cholesky_solver
{
public:
	static const bool is_serial = true;
	using stypes = cdr2::stypes;
	using stencil_op = typename stypes::template stencil_op<nine_pt>;
	using grid_func = typename stypes::grid_func;
	cholesky_solver(stencil_op & sop, std::shared_ptr<config> conf);
	template<class T>
	cholesky_solver(stencil_op & sop, std::shared_ptr<config> conf, T & par_sman) : conf(conf)
	{
		setup(sop);
	}
	~cholesky_solver();
	void cycle(grid_func & x, const grid_func & b);
	config & get_config() { return *conf; }

protected:
	kman_ptr kman;
	grid_func ABD;
	real_t *bbd;
	std::shared_ptr<config> conf;
	void setup(stencil_op & sop);
};

}}

#endif
