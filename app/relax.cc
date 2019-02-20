#include <mpi.h>
#include <random>
#include <limits>

#include <cedar/3d/util/topo.h>
#include <cedar/3d/mpi/gallery.h>
#include <cedar/3d/mpi/kernel_manager.h>
#include <cedar/cycle/types.h>


static void set_random(cedar::cdr3::mpi::grid_func & x)
{
	using namespace cedar;

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::mt19937 gen;
	gen.seed(rank);
	std::uniform_real_distribution<real_t> dis;

	for (auto k : x.range(2)) {
		for (auto j : x.range(1)) {
			for (auto i : x.range(0)) {
				x(i,j,k) = dis(gen);
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

	config conf("openmp-relax.json");
	log::init(conf);

	auto kreg = mpi::build_kernel_manager(conf);
	auto topo = util::create_topo(conf);

	auto so = mpi::gallery::fe(topo);

	auto b = mpi::grid_func::ones(topo);
	auto x = mpi::grid_func::zeros(topo);

	set_random(x);

	// setup halo
	{
		auto & halo_service = kreg->services().get<mpi::halo_exchange>();
		std::vector<topo_ptr> topos{{so.grid_ptr()}};
		halo_service.setup(topos);
		halo_service.run(so);
	}

	relax_stencil sor(topo->nlocal(0) - 2,
	                  topo->nlocal(1) - 2,
	                  topo->nlocal(2) - 2);
	kreg->setup<mpi::point_relax>(so, sor);

	int nsweeps = conf.get<int>("nsweeps");
	for (auto i : range(nsweeps)) {
		(void)i;
		kreg->run<mpi::point_relax>(so, x, b, sor, cycle::Dir::DOWN);
	}

	auto res = mpi::grid_func::zeros(topo);
	kreg->run<mpi::residual>(so, x, b, res);

	real_t nrm = res.lp_norm<2>();
	std::cout.precision(std::numeric_limits<real_t>::max_digits10);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
		std::cout << "residual norm: " << nrm << std::endl;

	MPI_Finalize();
	return 0;
}
