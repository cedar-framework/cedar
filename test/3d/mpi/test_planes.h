#ifndef CEDAR_TEST_3D_MPI_PLANES_H
#define CEDAR_TEST_3D_MPI_PLANES_H

#include <gtest/gtest.h>
#include <cedar/3d/gallery.h>
#include <cedar/3d/mpi/gallery.h>
#include <cedar/3d/util/topo.h>
#include <cedar/3d/types.h>
#include <cedar/3d/mpi/types.h>
#include <cedar/3d/kernel_manager.h>
#include <cedar/3d/mpi/kernel_manager.h>

void set_rhs(cedar::cdr3::grid_func & b);
void set_rhs(cedar::cdr3::mpi::grid_func & b);
void set_sol(cedar::cdr3::grid_func & x);
void set_sol(cedar::cdr3::mpi::grid_func & x);


template<class sten>
cedar::cdr3::stencil_op<sten> create_op_ser(cedar::len_t nx, cedar::len_t ny, cedar::len_t nz);

template<class sten>
cedar::cdr3::mpi::stencil_op<sten> create_op_mpi(cedar::topo_ptr grid);


template<cedar::relax_dir rdir, class sten>
void plane_test(cedar::len_t nx, cedar::len_t ny, cedar::len_t nz)
{
	using namespace cedar;
	using namespace cedar::cdr3;

	config::reader conf("test-planes-mpi.json");
	log::init(conf);

	int nsweeps = 3;

	auto kman_ser = build_kernel_manager(conf);
	auto kman_mpi = mpi::build_kernel_manager(conf);

	auto grid = util::create_topo_global(MPI_COMM_WORLD, nx, ny, nz);
	auto so_ser = create_op_ser<sten>(nx, ny, nz);
	grid_func b_ser(nx, ny, nz), x_ser(nx, ny, nz);
	set_rhs(b_ser);
	set_sol(x_ser);

	auto so_mpi = create_op_mpi<sten>(grid);
	mpi::grid_func b_mpi(grid), x_mpi(grid);
	set_rhs(b_mpi);
	set_sol(x_mpi);

	// setup halo
	{
		std::vector<topo_ptr> topos{{so_mpi.grid_ptr()}};
		kman_mpi->setup<mpi::halo_exchange>(topos);
		kman_mpi->run<mpi::halo_exchange>(so_mpi);
		kman_mpi->run<mpi::halo_exchange>(x_mpi);
		kman_mpi->run<mpi::halo_exchange>(b_mpi);
	}

	kman_ser->setup<plane_relax<rdir>>(so_ser);
	kman_mpi->setup<mpi::plane_relax<rdir>>(so_mpi);

	for (auto i : range(nsweeps)) {
		(void)i;
		kman_ser->run<plane_relax<rdir>>(so_ser, x_ser, b_ser, cycle::Dir::DOWN);
		kman_ser->run<plane_relax<rdir>>(so_ser, x_ser, b_ser, cycle::Dir::UP);
		kman_mpi->run<mpi::plane_relax<rdir>>(so_mpi, x_mpi, b_mpi, cycle::Dir::DOWN);
		kman_mpi->run<mpi::plane_relax<rdir>>(so_mpi, x_mpi, b_mpi, cycle::Dir::UP);
	}

	real_t tol = 1e-7;

	for (auto k : x_mpi.range(2)) {
		for (auto j : x_mpi.range(1)) {
			for (auto i : x_mpi.range(0)) {
				len_t is = grid->is(0) + i - 1;
				len_t js = grid->is(1) + j - 1;
				len_t ks = grid->is(2) + k - 1;

				auto diff = x_mpi(i,j,k) - x_ser(is,js,ks);
				ASSERT_LT(std::abs(diff), tol);
			}
		}
	}
}


#endif
