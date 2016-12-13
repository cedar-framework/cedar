#include <math.h>
#include <memory>
#include <iostream>

#include <boxmg/types.h>
#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/solver.h>

extern "C" {
	using namespace boxmg;
	void putf(real_t *so, real_t *qf,
	          len_t ii, len_t jj, len_t kk,
	          real_t hx, real_t hy, real_t hz);
	void put_sol(real_t *q,
	             len_t ii, len_t jj, len_t kk,
	             real_t hx, real_t hy, real_t hz);
}


static void set_problem(boxmg::bmg3::stencil_op & so, boxmg::bmg3::grid_func & b)
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	const double pi = M_PI;

	auto nx = b.shape(0);
	auto ny = b.shape(1);
	auto nz = b.shape(2);

	real_t hx = 1.0/(nx+1);
	real_t hy = 1.0/(ny+1);
	real_t hz = 1.0/(nz+1);

	auto & sten = so.stencil();
	sten.five_pt() = true;
	sten.set(0);
	b.set(0);

	putf(so.data(), b.data(),
	     b.len(0), b.len(1), b.len(2),
	     hx, hy, hz);
}


static void set_solution(boxmg::bmg3::grid_func & x)
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	auto nx = x.shape(0);
	auto ny = x.shape(1);
	auto nz = x.shape(2);

	real_t hx = 1.0/(nx+1);
	real_t hy = 1.0/(ny+1);
	real_t hz = 1.0/(nz+1);

	put_sol(x.data(),
	        x.len(0), x.len(1), x.len(2),
	        hx, hy, hz);
}


int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	log::status << "Beginning test" << std::endl;

	config::reader conf;

	auto ndofs = conf.getvec<len_t>("grid.n");
	auto nx = ndofs[0];
	auto ny = ndofs[1];
	auto nz = ndofs[2];

	auto so = stencil_op(nx, ny, nz);
	grid_func b(nx, ny, nz);
	set_problem(so, b);

	solver bmg(std::move(so));

	auto sol = bmg.solve(b);

	grid_func exact_sol(nx, ny, nz);

	set_solution(exact_sol);

	auto diff = exact_sol - sol;

	log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	// {
	// 	std::ofstream sten_file;
	// 	sten_file.open("stencil", std::ios::out | std::ios::trunc | std::ios::binary);
	// 	sten_file << bmg.level(-1).A;
	// 	sten_file.close();
	// }

	log::status << "Finished test" << std::endl;

	return 0;
}
