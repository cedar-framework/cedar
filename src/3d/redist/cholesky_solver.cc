#include <cedar/3d/redist/cholesky_solver.h>
#include <cedar/3d/mpi/grid_func.h>

namespace cedar { namespace cdr3 {

cholesky_solver::cholesky_solver(stencil_op & sop,
                                 std::shared_ptr<config::reader> conf) : conf(conf)
{
	this->kman = build_kernel_manager(*conf);
	auto nxc = sop.shape(0);
	auto nyc = sop.shape(1);
	auto nzc = sop.shape(2);
	this->ABD = mpi::grid_func(nxc*(nyc+1)+2, nxc*nyc*nzc, 0);
	this->bbd = new real_t[this->ABD.len(1)];

	kman->setup<solve_cg>(sop, ABD);
}


void cholesky_solver::cycle(grid_func & x, const grid_func & b)
{
	kman->run<solve_cg>(x, b, ABD, bbd);
}

}}
