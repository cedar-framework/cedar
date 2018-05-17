#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/mpi/gallery.h>

#include <cedar/util/time_log.h>

using namespace cedar;
using namespace cedar::cdr2;

static mpi::stencil_op<five_pt> create_op(topo_ptr grid, std::array<bool, 3> periodic)
{
	auto so = mpi::stencil_op<five_pt>(grid);
	auto & topo = so.grid();

	so.set(0);

	real_t nlx = topo.nlocal(0) - 2;
	real_t nly = topo.nlocal(1) - 2;
	real_t ngx = topo.nglobal(0) - 2;
	real_t ngy = topo.nglobal(1) - 2;

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);

	auto nx = topo.nglobal(0) - 2;
	auto ny = topo.nglobal(1) - 2;

	if (periodic[0]) nx--;
	if (periodic[1]) ny--;

	real_t hx = 1.0 / (nx + 1);
	real_t hy = 1.0 / (ny + 1);

	real_t nlx_g = nlx + 2;
	real_t nly_g = nly + 2;

	real_t i1 = nlx + 1;
	real_t j1 = nly + 1;
	real_t i2 = nlx;
	real_t j2 = nly;

	real_t igf = igs + i2 - 1;
	real_t jgf = jgs + j2 - 1;

	real_t xh = hy / hx;
	real_t yh = hx / hy;

	real_t ibeg = 1;
	real_t jbeg = 1;

	if (igs == 1 and not periodic[0])
		ibeg++;
	if (jgs == 1 and not periodic[1])
		jbeg++;

	real_t iend = nlx_g;
	real_t jend = nly_g;

	if (igf == ngx and not periodic[0])
		iend--;
	if (jgf == ngy and not periodic[1])
		jend--;

	for (auto j : range<len_t>(jbeg, jend)) {
		for (auto i : range<len_t>(1, iend)) {
			so(i,j,five_pt::s) = yh;
		}
	}

	for (auto j : range<len_t>(1, jend)) {
		for (auto i : range<len_t>(ibeg, iend)) {
			so(i,j,five_pt::w) = xh;
		}
	}

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			so(i,j,five_pt::c) = 2*xh + 2*yh;
		}
	}

	return so;
}


static void set_problem(mpi::grid_func & b, std::array<bool, 3> periodic)
{
	const double pi = M_PI;

	auto rhs = [pi](real_t x, real_t y) {
		return 8*(pi*pi)*sin(2*pi*x)*sin(2*pi*y);
	};

	auto & topo = b.grid();

	b.set(0);

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);

	auto ngx = topo.nglobal(0) - 2;
	auto ngy = topo.nglobal(1) - 2;

	if (periodic[0]) ngx--;
	if (periodic[1]) ngy--;

	real_t hx = 1.0 / (ngx + 1);
	real_t hy = 1.0 / (ngy + 1);

	real_t h2 = hx*hy;

	real_t nlx = topo.nlocal(0) - 2;
	real_t nly = topo.nlocal(1) - 2;

	real_t i1 = nlx + 1;
	real_t j1 = nly + 1;

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			len_t is = igs + i;
			len_t js = jgs + j;

			real_t x = (is-1)*hx;
			real_t y = (js-1)*hy;

			b(i,j) = rhs(x, y) * h2;

		}
	}
}


static void set_solution(mpi::grid_func & q, std::array<bool, 3> periodic)
{
	const double pi = M_PI;

	auto sol = [pi](real_t x, real_t y) {
		return sin(2*pi*x)*sin(2*pi*y);
	};

	auto & topo = q.grid();

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);

	auto ngx = topo.nglobal(0) - 2;
	auto ngy = topo.nglobal(1) - 2;

	if (periodic[0]) ngx--;
	if (periodic[1]) ngy--;

	real_t hx = 1.0 / (ngx + 1);
	real_t hy = 1.0 / (ngy + 1);

	for (auto j : q.range(1)) {
		for (auto i : q.range(0)) {
			len_t is = igs + i;
			len_t js = jgs + j;

			real_t x = (is-1)*hx;
			real_t y = (js-1)*hy;

			q(i,j) = sol(x,y);
		}
	}
}


int main(int argc, char *argv[])
{
	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	timer_init(MPI_COMM_WORLD);

	auto conf = std::make_shared<config>();
	auto params = build_kernel_params(*conf);
	log::init(*conf);
	log::status << "Beginning test" << std::endl;
	auto grid = util::create_topo(*conf);

	mpi::grid_func b(grid);
	auto so = create_op(grid, params->periodic);

	set_problem(b, params->periodic);

	mpi::solver<five_pt> bmg(so, conf);

	// {
	// 	std::string suffix(std::to_string(grid->coord(0)) + "." +
	// 	                   std::to_string(grid->coord(1)));
	// 	std::ofstream ffile("output/stencil-" + suffix);
	// 	std::ofstream rfile("output/restrict-" + suffix);
	// 	std::ofstream cfile("output/coarse-" + suffix);
	// 	ffile << bmg.level(-1).A;
	// 	rfile << bmg.level(-1).P;
	// 	cfile << bmg.level(0).A;
	// 	rfile.close();
	// 	ffile.close();
	// 	cfile.close();
	// }

	auto sol = bmg.solve(b);

	mpi::grid_func exact_sol(sol.grid_ptr());
	set_solution(exact_sol, params->periodic);

	mpi::grid_func diff = exact_sol - sol;

	log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
