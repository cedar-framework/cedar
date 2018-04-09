#include <gtest/gtest.h>

#include <cedar/3d/util/topo.h>
#include <cedar/3d/mpi/gallery.h>
#include <cedar/3d/gallery.h>
#include <cedar/3d/kernel_manager.h>
#include <cedar/3d/mpi/kernel_manager.h>
#include <cedar/cycle/types.h>


TEST(MPIRelax3, Point7) {
	using namespace cedar;
	using namespace cedar::cdr3;

	config::reader conf("");
	log::init(conf);
	auto kreg_ser = build_kernel_manager(conf);
	auto kreg_mpi = mpi::build_kernel_manager(conf);

	int nsweeps = 5;

	auto nx = 50;
	auto ny = nx;
	auto nz = nx;

	auto grid = util::create_topo_global(MPI_COMM_WORLD, nx, ny, nz);
	auto so_mpi = mpi::gallery::poisson(grid);
	auto b_mpi = mpi::grid_func::zeros(grid);
	auto x_mpi = mpi::grid_func::ones(grid);

	auto so_ser = gallery::poisson(nx, ny, nz);
	auto b_ser = grid_func::zeros(nx, ny, nz);
	auto x_ser = grid_func::ones(nx, ny, nz);

	// setup halo
	{
		std::vector<topo_ptr> topos{{so_mpi.grid_ptr()}};
		kreg_mpi->setup<mpi::halo_exchange>(topos);
		kreg_mpi->run<mpi::halo_exchange>(so_mpi);
	}

	relax_stencil sor_mpi(nx, ny, nz);
	relax_stencil sor_ser(nx, ny, nz);

	kreg_ser->setup<point_relax>(so_ser, sor_ser);
	kreg_mpi->setup<mpi::point_relax>(so_mpi, sor_mpi);

	for (auto i : range(nsweeps)) {
		(void)i;
		kreg_ser->run<point_relax>(so_ser, x_ser, b_ser, sor_ser, cycle::Dir::DOWN);
		kreg_ser->run<point_relax>(so_ser, x_ser, b_ser, sor_ser, cycle::Dir::UP);
		kreg_mpi->run<mpi::point_relax>(so_mpi, x_mpi, b_mpi, sor_mpi, cycle::Dir::DOWN);
		kreg_mpi->run<mpi::point_relax>(so_mpi, x_mpi, b_mpi, sor_mpi, cycle::Dir::UP);
	}

	auto ndiff = x_mpi.lp_norm<2>() - x_ser.lp_norm<2>();

	real_t tol = 1e-10;
	ASSERT_LT(std::abs(ndiff), tol);
}
