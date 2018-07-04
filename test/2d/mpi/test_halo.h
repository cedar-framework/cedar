#ifndef CEDAR_TEST_2D_MPI_HALO_H
#define CEDAR_TEST_2D_MPI_HALO_H

#include <gtest/gtest.h>
#include <mpi.h>

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/mpi/kernel_manager.h>
#include <cedar/2d/util/topo.h>
#include <cedar/kernel_params.h>

namespace cedar { namespace cdr2 {


static void fill_gfunc(mpi::grid_func & b)
{
	auto & topo = b.grid();
	for (auto j : b.range(1)) {
		for (auto i : b.range(0)) {
			len_t ci = (topo.is(0)-1) + (i-1);
			len_t cj = (topo.is(1)-1) + (j-1);
			b(i,j) = cj * topo.nglobal(0) + ci;
		}
	}
}


template<class sten>
static void fill_stencil(mpi::stencil_op<sten> & so)
{
	auto & topo = so.grid();
	for (auto j : so.range(1)) {
		for (auto i : so.range(0)) {
			for (int k = 0; k < stencil_ndirs<sten>::value; k++) {
				len_t ci = (topo.is(0)-1) + (i-1);
				len_t cj = (topo.is(1)-1) + (j-1);
				so(i,j,static_cast<sten>(k)) = cj * topo.nglobal(0) * stencil_ndirs<sten>::value
					+ ci * stencil_ndirs<sten>::value + k;
			}
		}
	}
}


static bool set_coords_periodic(int & ci, int & cj,
                                len_t i, len_t j,
                                len_t lenx, len_t leny,
                                grid_topo & topo, kernel_params & params)
{
	if ((i == 0) and (topo.coord(0) == 0)) {
		if (not params.periodic[0])
			return false;
		ci = topo.nglobal(0) - 3;
	} else if ((i == (lenx - 1)) and (topo.coord(0) == (topo.nproc(0) - 1))) {
		if (not params.periodic[0])
			return false;
		ci = 0;
	} else
		ci = (topo.is(0)-1) + (i-1);

	if ((j == 0) and (topo.coord(1) == 0)) {
		if (not params.periodic[1])
			return false;
		cj = topo.nglobal(1) - 3;
	} else if ((j == (leny - 1)) and (topo.coord(1) == (topo.nproc(1) - 1))) {
		if (not params.periodic[1])
			return false;
		cj = 0;
	} else
		cj = (topo.is(1)-1) + (j-1);

	return true;
}


static void test_gfunc(kernel_params & params, mpi::grid_func & b)
{
	auto & topo = b.grid();
	for (auto j : b.grange(1)) {
		for (auto i : b.grange(0)) {
			int ci, cj;
			if (not set_coords_periodic(ci, cj, i, j, b.len(0), b.len(1), topo, params))
				continue;

			/* if (b(i,j) != cj * topo.nglobal(0) + ci) */
			/* 	printf("[%d,%d] %d %d\n", topo.coord(0), topo.coord(1), i, j); */
			ASSERT_EQ(b(i,j), cj * topo.nglobal(0) + ci);
		}
	}

}


template<class sten>
static void test_stencil(kernel_params & params, mpi::stencil_op<sten> & so)
{
	auto & topo = so.grid();
	for (auto j : so.grange(1)) {
		for (auto i : so.grange(0)) {
			for (int k = 0; k < stencil_ndirs<sten>::value; k++) {
				int ci, cj;
				if (not set_coords_periodic(ci, cj, i, j, so.len(0), so.len(1), topo, params))
					continue;

				ASSERT_EQ(so(i,j,static_cast<sten>(k)),
				          cj * topo.nglobal(0) * stencil_ndirs<sten>::value
				          + ci * stencil_ndirs<sten>::value + k);
			}
		}
	}
}


static void set_conf(config & conf, const std::array<int, 2> & proc, int per_mask)
{
	conf.set("grid.local", true);
	std::vector<int> gsize{{20, 30}};
	std::vector<int> pgrid{{proc[0], proc[1]}};
	conf.setvec("grid.n", gsize);
	conf.setvec("grid.np", pgrid);
	std::vector<bool> periodic{{false, false}};
	for (int i = 0; i < 3; i++) {
		if ((1<<i) & per_mask) periodic[i] = true;
	}
	conf.setvec("grid.periodic", periodic);
}


void run_test(const std::string & halo_name, MPI_Comm comm, std::array<int, 2> & proc, int per_mask)
{
	auto conf = std::make_shared<config>("config.json");
	log::set_comm(comm);
	set_conf(*conf, proc, per_mask);
	auto params = build_kernel_params(*conf);
	log::init(*conf);

	auto grid = util::create_topo(*conf, comm);
	mpi::grid_func b(grid);
	mpi::stencil_op<five_pt> so(grid);

	fill_gfunc(b);
	fill_stencil(so);
	mpi::solver<five_pt> slv(so, conf);
	auto kman = slv.get_kernels();
	auto & halo_service = kman->services().get<mpi::halo_exchange>();
	halo_service.run(b);
	halo_service.run(so);
	test_gfunc(*params, b);
	test_stencil(*params, so);

	for (std::size_t lvl = 1; lvl < slv.nlevels(); lvl++) {
		auto & level = slv.levels.get(lvl);
		fill_gfunc(level.b);
		fill_stencil(level.A);
		halo_service.run(level.b);
		halo_service.run(level.A);
		test_gfunc(*params, level.b);
		test_stencil(*params, level.A);
	}
}


void test_driver(const std::string & halo_name, std::vector<std::array<int,2>> & procs)
{
	for (auto & proc : procs) {
		int np = proc[0]*proc[1];
		int rank;
		MPI_Comm comm;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_split(MPI_COMM_WORLD, rank < np, rank, &comm);
		if (rank < np) {
			for (int per_mask = 0; per_mask < 4; per_mask++) {
				run_test(halo_name, comm, proc, per_mask);
			}
		}
	}
}

}}

#endif
