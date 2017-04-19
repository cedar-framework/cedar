#include <mpi.h>
#include <iostream>
#include <memory>
#include <array>
#include <math.h>

#include <cedar/types.h>
#include <cedar/util/grid.h>
#include <cedar/2d/grid_func.h>
#include <cedar/2d/stencil_op.h>
#include <cedar/2d/solver.h>
#include <cedar/2d/gallery.h>


static cedar::cdr2::stencil_op create_op(cedar::len_t nx, cedar::len_t ny, std::array<bool, 3> periodic)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	auto so = stencil_op(nx, ny);
	auto & sten = so.stencil();
	sten.five_pt() = true;

	sten.set(0);

	if (periodic[0]) nx--;
	if (periodic[1]) ny--;
	real_t hx = 1.0/(nx+1);
	real_t hy = 1.0/(ny+1);
	real_t xh = hy/hx;
	real_t yh = hx/hy;
	len_t l = sten.shape(0);
	len_t m = sten.shape(1);
	len_t i1 = sten.shape(0)+1;
	len_t j1 = sten.shape(1)+1;
	len_t ibeg = 2;
	len_t jbeg = 2;

	if (periodic[0]) ibeg--;
	if (periodic[1]) jbeg--;

	auto & o = sten;

	for (auto j : range<len_t>(jbeg, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			o(i,j,dir::S) = 1.0 * yh;
		}
	}

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(ibeg, i1)) {
			o(i,j,dir::W) = 1.0 * xh;
		}
	}

	for (auto j : sten.range(1)) {
		for (auto i : sten.range(0)) {
			o(i,j,dir::C) = 2*xh + 2*yh;
		}
	}

	if (periodic[0]) {
		for (auto j : sten.range(1)) {
			o(ibeg-1,j,dir::C) = o(l,j,dir::C);
			o(ibeg-1,j,dir::W) = o(l,j,dir::W);
			o(ibeg-1,j,dir::S) = o(l,j,dir::S);

			o(l+1,j,dir::C) = o(ibeg,j,dir::C);
			o(l+1,j,dir::W) = o(ibeg,j,dir::W);
			o(l+1,j,dir::S) = o(ibeg,j,dir::S);
		}
	}

	if (periodic[1]) {
		for (auto i : sten.range(0)) {
			o(i,jbeg-1,dir::C) = o(i,m,dir::C);
			o(i,jbeg-1,dir::W) = o(i,m,dir::W);
			o(i,jbeg-1,dir::S) = o(i,m,dir::S);

			o(i,m+1,dir::C) = o(i,jbeg,dir::C);
			o(i,m+1,dir::W) = o(i,jbeg,dir::W);
			o(i,m+1,dir::S) = o(i,jbeg,dir::S);
		}
	}

	return so;
}


static void set_problem(cedar::cdr2::grid_func & b, std::array<bool, 3> periodic)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	const double pi = M_PI;

	auto rhs = [pi](real_t x, real_t y) {
		return 8*(pi*pi)*sin(2*pi*x)*sin(2*pi*y);
	};

	b.set(0);

	len_t nx = b.len(0) - 2;
	len_t ny = b.len(1) - 2;
	if (periodic[0]) nx--;
	if (periodic[1]) ny--;

	real_t hx = 1.0/(nx + 1);
	real_t hy = 1.0/(ny + 1);
	real_t h2 = hx*hy;

	for (auto j : b.range(1)) {
		for (auto i : b.range(0)) {
			real_t x = i*hx;
			real_t y = j*hy;

			b(i,j) = rhs(x,y) * h2;
		}
	}

	if (periodic[0]) {
		for (auto j : b.grange(1)) {
			b(0,j) = b(b.shape(0),j);
			b(b.shape(0)+1,j) = b(1,j);
		}
	}

	if (periodic[1]) {
		for (auto i : b.grange(0)) {
			b(i,0) = b(i,b.shape(1));
			b(i,b.shape(1)+1) = b(i,1);
		}
	}
}


static void set_solution(cedar::cdr2::grid_func & q, std::array<bool, 3> periodic)
{
	using namespace cedar;

	const double pi = M_PI;

	auto sol = [pi](real_t x, real_t y) {
		return sin(2*pi*x)*sin(2*pi*y);
	};

	len_t nx = q.len(0) - 2;
	len_t ny = q.len(1) - 2;
	if (periodic[0]) nx--;
	if (periodic[1]) ny--;

	real_t hx = 1.0/(nx + 1);
	real_t hy = 1.0/(ny + 1);
	for (auto j : q.range(1)) {
		for (auto i : q.range(0)) {
			real_t x = i*hx;
			real_t y = j*hy;

			q(i,j) = sol(x, y);
		}
	}
}


int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr2;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	auto conf = std::make_shared<config::reader>();
	auto params = build_kernel_params(*conf);
	auto ndofs = conf->getvec<len_t>("grid.n");
	auto nx = ndofs[0];
	auto ny = ndofs[1];

	auto so = create_op(nx, ny, params->periodic);
	grid_func b(nx, ny);

	set_problem(b, params->periodic);

	solver bmg(std::move(so), conf);
	{
		std::ofstream ffile("fine.txt");
		std::ofstream cfile("coarse.txt");
		std::ofstream rfile("restrict.txt");
		ffile << bmg.level(-1).A;
		cfile << bmg.level(-2).A;
		rfile << bmg.level(-1).P;
		ffile.close();
		cfile.close();
		rfile.close();
	}


	auto sol = bmg.solve(b);

	grid_func exact_sol(sol.shape(0), sol.shape(1));

	// std::ofstream dfile;
	// dfile.open("sol.txt", std::ios::out | std::ios::trunc | std::ios::binary);

	set_solution(exact_sol, params->periodic);

	auto diff = exact_sol - sol;
	// dfile << diff << std::endl;

	// dfile.close();

	log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
