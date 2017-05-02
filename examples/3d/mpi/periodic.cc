#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <cedar/types.h>
#include <cedar/3d/mpi/grid_func.h>
#include <cedar/3d/util/topo.h>
#include <cedar/3d/mpi/solver.h>
#include <cedar/3d/mpi/gallery.h>

#include <cedar/util/time_log.h>


static cedar::cdr3::mpi::stencil_op create_op(cedar::topo_ptr grid, std::array<bool, 3> periodic)
{
	using namespace cedar;
	using namespace cedar::cdr3;

	auto so = mpi::stencil_op(grid);
	auto & o = so.stencil();
	o.five_pt() = true;

	auto & topo = so.grid();

	o.set(0);

	real_t nlx = topo.nlocal(0) - 2;
	real_t nly = topo.nlocal(1) - 2;
	real_t nlz = topo.nlocal(2) - 2;
	real_t ngx = topo.nglobal(0) - 2;
	real_t ngy = topo.nglobal(1) - 2;
	real_t ngz = topo.nglobal(2) - 2;

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);
	real_t kgs = topo.is(2);

	auto nx = topo.nglobal(0) - 2;
	auto ny = topo.nglobal(1) - 2;
	auto nz = topo.nglobal(2) - 2;

	if (periodic[0]) nx--;
	if (periodic[1]) ny--;
	if (periodic[2]) nz--;

	real_t hx = 1.0 / (nx + 1);
	real_t hy = 1.0 / (ny + 1);
	real_t hz = 1.0 / (nz + 1);

	real_t nlx_g = nlx + 2;
	real_t nly_g = nly + 2;
	real_t nlz_g = nlz + 2;

	real_t i1 = nlx + 1;
	real_t j1 = nly + 1;
	real_t k1 = nlz + 1;
	real_t i2 = nlx;
	real_t j2 = nly;
	real_t k2 = nlz;

	real_t igf = igs + i2 - 1;
	real_t jgf = jgs + j2 - 1;
	real_t kgf = kgs + k2 - 1;

	real_t xh=hy*hz/hx;
	real_t yh=hx*hz/hy;
	real_t zh=hx*hy/hz;

	real_t ibeg = 1;
	real_t jbeg = 1;
	real_t kbeg = 1;

	if (igs == 1 and not periodic[0])
		ibeg++;
	if (jgs == 1 and not periodic[1])
		jbeg++;
	if (kgs == 1 and not periodic[2])
		kbeg++;

	real_t iend = nlx_g;
	real_t jend = nly_g;
	real_t kend = nlz_g;

	if (igf == ngx and not periodic[0])
		iend--;
	if (jgf == ngy and not periodic[1])
		jend--;
	if (kgf == ngz and not periodic[2])
		kend--;

	for (auto k : range<len_t>(1, kend)) {
		for (auto j : range<len_t>(jbeg, jend)) {
			for (auto i : range<len_t>(1, iend)) {
				o(i,j,k,dir::PS) = 1.0 * yh;
			}
		}
	}

	for (auto k : range<len_t>(1, kend)) {
		for (auto j : range<len_t>(1, jend)) {
			for (auto i : range<len_t>(ibeg, iend)) {
				o(i,j,k,dir::PW) = 1.0 * xh;
			}
		}
	}

	for (auto k : range<len_t>(kbeg, kend)) {
		for (auto j : range<len_t>(1, jend)) {
			for (auto i : range<len_t>(1, iend)) {
				o(i,j,k,dir::B) = 1.0 * zh;
			}
		}
	}

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				o(i,j,k,dir::P) = 2*xh + 2*yh + 2*zh;
			}
		}
	}

	return so;
}


static void set_problem(cedar::cdr3::mpi::grid_func & b, std::array<bool, 3> periodic)
{
	using namespace cedar;
	using namespace cedar::cdr3;

	const double pi = M_PI;

	auto rhs = [pi](real_t x, real_t y, real_t z) {
		return 12*(pi*pi)*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z);
	};

	auto & topo = b.grid();

	b.set(0);

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);
	real_t kgs = topo.is(2);

	auto ngx = topo.nglobal(0) - 2;
	auto ngy = topo.nglobal(1) - 2;
	auto ngz = topo.nglobal(2) - 2;

	if (periodic[0]) ngx--;
	if (periodic[1]) ngy--;
	if (periodic[2]) ngz--;

	real_t hx = 1.0 / (ngx + 1);
	real_t hy = 1.0 / (ngy + 1);
	real_t hz = 1.0 / (ngz + 1);

	real_t h2 = hx*hy*hz;

	real_t nlx = topo.nlocal(0) - 2;
	real_t nly = topo.nlocal(1) - 2;
	real_t nlz = topo.nlocal(2) - 2;

	real_t i1 = nlx + 1;
	real_t j1 = nly + 1;
	real_t k1 = nlz + 1;

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				len_t is = igs + i;
				len_t js = jgs + j;
				len_t ks = kgs + k;

				real_t x = (is-1)*hx;
				real_t y = (js-1)*hy;
				real_t z = (ks-1)*hz;

				b(i,j,k) = rhs(x, y, z) * h2;

			}
		}
	}
}


int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr3;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	timer_init(MPI_COMM_WORLD);

	auto conf = std::make_shared<config::reader>();
	auto params = build_kernel_params(*conf);
	log::init(*conf);
	log::status << "Beginning test" << std::endl;
	auto grid = util::create_topo(*conf);

	mpi::grid_func b(grid);
	auto so = create_op(grid, params->periodic);

	set_problem(b, params->periodic);

	mpi::solver bmg(std::move(so), conf);

	{
		std::string suffix(std::to_string(grid->coord(0)) + "." +
		                   std::to_string(grid->coord(1)) + "." +
		                   std::to_string(grid->coord(2)));
		std::ofstream ffile("output/fine-" + suffix);
		std::ofstream rfile("output/restrict-" + suffix);
		std::ofstream cfile("output/coarse-" + suffix);
		ffile << bmg.level(-1).A;
		rfile << bmg.level(-1).P;
		cfile << bmg.level(-2).A;
		rfile.close();
		ffile.close();
		cfile.close();
	}

	// auto sol = bmg.solve(b);

	// mpi::grid_func exact_sol(sol.grid_ptr());
	// set_solution(exact_sol, params->periodic);

	// mpi::grid_func diff = exact_sol - sol;

	// log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
