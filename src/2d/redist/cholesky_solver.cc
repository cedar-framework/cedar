#include <cedar/2d/redist/cholesky_solver.h>

namespace cedar { namespace cdr2 {

cholesky_solver::cholesky_solver(stencil_op & sop,
                                 std::shared_ptr<config> conf) : conf(conf)
{
	setup(sop);
}


void cholesky_solver::setup(stencil_op & sop)
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

        std::cerr << " :: Setting up Cholesky solver." << std::endl;

        sop.ensure_cpu();
	kman->setup<solve_cg>(sop, ABD);
}


cholesky_solver::~cholesky_solver()
{
	if (bbd)
		delete[] bbd;
}


void cholesky_solver::cycle(grid_func & x, const grid_func & b)
{
    std::cerr << " :: Calling Cholesky solve" << std::endl;
    auto& bd = const_cast<grid_func&>(b);
    x.ensure_cpu();
    bd.ensure_cpu();
    kman->run<solve_cg>(x, b, ABD, bbd);
    x.mark_cpu_dirty();
}

}}
