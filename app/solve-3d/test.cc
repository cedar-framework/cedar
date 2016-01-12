#include <math.h>
#include <memory>
#include <iostream>

#include <boxmg/types.h>
#include <boxmg/util/grid.h>
#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/solver.h>

extern "C" {
	using namespace boxmg;
	void putf(real_t *so, real_t *qf,
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


int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	log::status << "Beginning test" << std::endl;

	auto nx = config::get<len_t>("grid.nx", 301);
	auto ny = config::get<len_t>("grid.ny", 301);
	auto nz = config::get<len_t>("grid.nz", 301);

	auto so = stencil_op(nx, ny, nz);
	grid_func b(nx, ny, nz);
	set_problem(so, b);

	log::status << nx << " " << so.stencil().shape(0) << " " << so.stencil().len(0) << std::endl;

	//solver bmg(std::move(so));

	// auto sol = bmg.solve(b);
	log::status << so.stencil()(1, 1, 1, dir::P) << std::endl;
	std::ofstream sten_file;
	sten_file.open("Stencil", std::ios::out | std::ios::trunc | std::ios::binary);
	sten_file << so;
	sten_file.close();

	log::status << "Finished test" << std::endl;

	return 0;
}
