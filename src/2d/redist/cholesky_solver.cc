#include <cedar/2d/redist/cholesky_solver.h>

namespace cedar { namespace cdr2 {

cholesky_solver::cholesky_solver(stencil_op & sop,
                                 std::shared_ptr<config::reader> conf) : conf(conf)
{
	this->kman = build_kernel_manager(*conf);
	auto params = build_kernel_params(*conf);
	auto nxc = sop.shape(0);
	auto nyc = sop.shape(1);

	len_t abd_len_0 = nxc + 2;
	if (params->periodic[0] or params->periodic[1])
		abd_len_0 = nxc*nyc;
	this->ABD = grid_func(abd_len_0, nxc*nyc, 0);
	this->bbd = new real_t[this->ABD.len(1)];

	kman->setup<solve_cg>(sop, ABD);
}


void cholesky_solver::cycle(grid_func & x, const grid_func & b)
{
	kman->run<solve_cg>(x, b, ABD, bbd);
}

}}
