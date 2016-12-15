#include <gtest/gtest.h>

#include <boxmg/3d/util/topo.h>
#include <boxmg/3d/mpi/gallery.h>
#include <boxmg/3d/gallery.h>
#include <boxmg/cycle/types.h>


TEST(MPIRelax3, Point7) {
	using namespace boxmg;
	using namespace boxmg::bmg3;

	int nsweeps = 1;

	auto nx = 5;
	auto ny = nx;
	auto nz = nx;

	auto grid = util::create_topo_global(MPI_COMM_WORLD, nx, ny, nz);
	auto so_mpi = mpi::gallery::poisson(grid);
	auto b_mpi = mpi::grid_func::zeros(grid);
	auto x_mpi = mpi::grid_func::ones(grid);

	auto so_ser = gallery::poisson(nx, ny, nz);
	auto b_ser = grid_func::zeros(nx, ny, nz);
	auto x_ser = grid_func::ones(nx, ny, nz);

	auto kreg_mpi = so_mpi.get_registry();
	auto kreg_ser = so_ser.get_registry();

	// setup halo
	{
		void *halo_ctx;
		kreg_mpi->halo_setup(so_mpi.grid(), &halo_ctx);
		so_mpi.halo_ctx = halo_ctx;
		kreg_mpi->halo_stencil_exchange(so_mpi);

		b_mpi.halo_ctx = halo_ctx;
		x_mpi.halo_ctx = halo_ctx;
	}

	relax_stencil sor_mpi(nx, ny, nz);
	relax_stencil sor_ser(nx, ny, nz);

	kreg_ser->setup_relax(so_ser, sor_ser);
	kreg_mpi->setup_relax(so_mpi, sor_mpi);

	for (auto i : range(nsweeps)) {
		(void)i;
		kreg_ser->relax(so_ser, x_ser, b_ser, sor_ser, cycle::Dir::DOWN);
		kreg_mpi->relax(so_mpi, x_mpi, b_mpi, sor_mpi, cycle::Dir::DOWN);
	}

	log::status << "SER: " << x_ser.lp_norm<2>() << std::endl;
	log::status << "MPI: " << x_mpi.lp_norm<2>() << std::endl;
}
