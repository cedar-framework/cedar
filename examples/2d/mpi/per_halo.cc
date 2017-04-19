#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/kernel/mpi/factory.h>
#include <cedar/2d/kernel/mpi/registry.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/util/mpi_grid.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/mpi/gallery.h>

#include <cedar/util/time_log.h>


static void draw(const cedar::cdr2::mpi::grid_func & b, std::ostream & os)
{
	for (auto j : b.grange(1)) {
		for (auto i : b.grange(0)) {
			if (b(i,j) < 0)
				os << '*';
			else
				os << b(i,j);
			os << " ";
		}
		os << '\n';
	}
}


static cedar::cdr2::mpi::stencil_op create_op(cedar::topo_ptr grid)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	auto so = mpi::stencil_op(grid);
	auto & o = so.stencil();
	o.five_pt() = true;

	auto & topo = so.grid();

	o.set(0);

	real_t nlx = topo.nlocal(0) - 2;
	real_t nly = topo.nlocal(1) - 2;
	real_t ngx = topo.nglobal(0) - 2;
	real_t ngy = topo.nglobal(1) - 2;

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);

	real_t hx = 1.0 / (topo.nglobal(0) - 1);
	real_t hy = 1.0 / (topo.nglobal(1) - 1);

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

	if (igs == 1)
		ibeg++;
	if (jgs == 1)
		jbeg++;

	real_t iend = nlx_g;
	real_t jend = nly_g;

	if (igf == ngx)
		iend--;
	if (jgf == ngy)
		jend--;

	for (auto j : range<len_t>(jbeg, jend)) {
		for (auto i : range<len_t>(1, iend)) {
			o(i,j,dir::S) = 1.0 * yh;
		}
	}

	for (auto j : range<len_t>(1, jend)) {
		for (auto i : range<len_t>(ibeg, iend)) {
			o(i,j,dir::W) = xh;
		}
	}

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			o(i,j,dir::C) = 2*xh + 2*yh;
		}
	}

	return so;
}

static void set_problem(cedar::cdr2::mpi::grid_func & b)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	const double pi = M_PI;

	auto rhs = [pi](real_t x, real_t y) {
		return 8*(pi*pi)*sin(2*pi*x)*sin(2*pi*y);
	};

	auto & topo = b.grid();

	b.set(0);

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);

	real_t hx = 1.0 / (topo.nglobal(0) - 1);
	real_t hy = 1.0 / (topo.nglobal(1) - 1);

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


static void set_solution(cedar::cdr2::mpi::grid_func & q)
{
	using namespace cedar;

	const double pi = M_PI;

	auto sol = [pi](real_t x, real_t y) {
		return sin(2*pi*x)*sin(2*pi*y);
	};

	auto & topo = q.grid();

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);

	real_t hx = 1.0 / (topo.nglobal(0) - 1);
	real_t hy = 1.0 / (topo.nglobal(1) - 1);

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
	using namespace cedar;
	using namespace cedar::cdr2;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	timer_init(MPI_COMM_WORLD);

	auto conf = std::make_shared<config::reader>();
	log::init(*conf);
	log::status << "Beginning test" << std::endl;
	auto grid = util::create_topo(*conf);


	mpi::grid_func b(grid);
	auto kreg = std::make_shared<cedar::cdr2::kernel::mpi::registry>();
	cedar::cdr2::kernel::mpi::factory::init(kreg, *conf);

	b.set(-1);
	void *halo_ctx;
	kreg->halo_setup(*grid, &halo_ctx);
	b.halo_ctx = halo_ctx;

	{
		for (auto j : b.range(1)) {
			for (auto i : b.range(0)) {
				b(i,j) = grid->coord(1)*grid->nproc(0) + grid->coord(0);
			}
		}

		std::ofstream ofile("output/before-" + std::to_string(grid->coord(0)) +
		                    "." + std::to_string(grid->coord(1)));
		draw(b, ofile);
		ofile.close();
	}

	kreg->halo_exchange(b);

	{
		std::ofstream ofile("output/after-" + std::to_string(grid->coord(0)) +
		                    "." + std::to_string(grid->coord(1)));
		draw(b, ofile);
		ofile.close();
	}

	// auto so = mpi::gallery::poisson(grid);
	// mpi::grid_func b(grid);

	// set_problem(b);

	// mpi::solver bmg(std::move(so));

	// MPI_Barrier(MPI_COMM_WORLD); // synchronize before timing solve
	// auto sol = bmg.solve(b);


	// mpi::grid_func exact_sol(sol.grid_ptr());
	// set_solution(exact_sol);

	// mpi::grid_func diff = exact_sol - sol;

	// log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	// timer_save("timings.json");

	log::status << "Finished Test" << std::endl;

	MPI_Finalize();
	return 0;
}
