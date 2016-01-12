#include <math.h>
#include <memory>
#include <iostream>

#include <boxmg/types.h>
#include <boxmg/util/grid.h>
#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/solver.h>


static void putf(boxmg::bmg3::stencil_op & so, boxmg::bmg3::grid_func & b,
          boxmg::len_t nx, boxmg::len_t ny, boxmg::len_t nz)
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	const double pi = M_PI;

	real_t hx = 1.0/(nx+1);
	real_t hy = 1.0/(ny+1);
	real_t hz = 1.0/(nz+1);

	real_t h2 = hx*hy*hz;
	real_t xh = hy*hz/hx;
	real_t yh = hx*hz/hy;
	real_t zh = hx*hy/hz;

	auto & sten = so.stencil();
	sten.five_pt() = true;
	sten.set(0);

	b.set(0);

	auto xs = linspace(0, 1, nx + 2);
	auto ys = linspace(0, 1, ny + 2);
	auto zs = linspace(0, 1, nz + 2);

	for (auto k : sten.range(2)) {
		for (auto j : range(static_cast<len_t>(2), sten.shape(1)+1)) {
			for (auto i : sten.range(0)) {
				sten(i, j, k, dir::PS) = 1.0 * yh;
			}
		}
	}

	for (auto k : sten.range(2)) {
		for (auto j : sten.range(1)) {
			for (auto i : range(static_cast<len_t>(2), sten.shape(0)+1)) {
				sten(i, j, k, dir::PW) = 1.0 * xh;
			}
		}
	}

	for (auto k : range(static_cast<len_t>(2), sten.shape(2)+1)) {
		for (auto j : sten.range(1)) {
			for (auto i : sten.range(0)) {
				sten(i, j, k, dir::B) = 1.0 * zh;
			}
		}
	}

	for (auto k : sten.range(2)) {
		for (auto j : sten.range(1)) {
			for (auto i: sten.range(0)) {
				sten(i, j, k, dir::P) = 2*yh + 2*zh + 2*xh;
				auto x = xs[i];
				auto y = ys[j];
				auto z = zs[k];
				b(i, j, k) = 12*(pi*pi)*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)*h2;
			}
		}
	}
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
	putf(so, b, nx, ny, nz);

	solver bmg(std::move(so));

	// auto sol = bmg.solve(b);
	std::ofstream sten_file;
	sten_file.open("Stencil", std::ios::out | std::ios::trunc | std::ios::binary);
	sten_file << bmg.level(-1).A;
	sten_file.close();

	log::status << "Finished test" << std::endl;

	return 0;
}
