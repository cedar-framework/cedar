#include <abt.h>

#include <cedar/2d/mpi/kernel_manager.h>
#include <cedar/2d/mpi/gallery.h>
#include <cedar/2d/util/topo.h>

#include <cedar/3d/mpi/plane_mempool.h>
#include <cedar/3d/mpi/plane_mpi.h>
#include <cedar/3d/mpi/plane_exchange.h>

extern "C" {
	void test_gather(cedar::real_t *tvals, int plane_len, cedar::real_t *recv, MPI_Fint comm, void *mp);
}

using namespace cedar;
using namespace cedar::cdr2;

struct ult_params
{
	mpi::stencil_op<five_pt> *A;
	mpi::grid_func *x;
	mpi::grid_func *b;
	relax_stencil *sor;
	mpi::grid_func *res;
	mpi::kman_ptr kman;
};

static void compute(void *args)
{
	ult_params *params = (ult_params*) args;

	auto kman = params->kman;
	auto & so = *params->A;
	auto & x = *params->x;
	auto & b = *params->b;
	auto & sor = *params->sor;
	auto &res = *params->res;

	kman->run<mpi::line_relax<relax_dir::x>>(so, x, b, sor, res, cycle::Dir::DOWN);
}

int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr2;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
	ABT_init(argc, argv);

	ABT_pool pool;
	ABT_xstream xstream;
	ABT_xstream_self(&xstream);
	ABT_xstream_get_main_pools(xstream, 1, &pool);

	config conf("config.json");
	log::init(conf);

	auto nplanes = conf.get<int>("nplanes");
	auto aggregate = conf.get<bool>("aggregate");

	auto grid = util::create_topo(conf);

	std::vector<mpi::kman_ptr> kmans;
	std::vector<mpi::grid_func> bs;
	std::vector<mpi::grid_func> xs;

	for (auto i : range<std::size_t>(nplanes)) {
		kmans.push_back(mpi::build_kernel_manager(conf));
		auto & sman = kmans.back()->services();
		if (i == 0) {
			sman.add<services::message_passing, cdr3::mpi::plane_setup_mpi>("plane-setup");
			sman.add<services::message_passing, cdr3::mpi::plane_mpi>("plane", nplanes);
			sman.add<services::mempool, cdr3::mpi::plane_mempool>("plane", nplanes);
			sman.add<mpi::halo_exchange, cdr3::mpi::plane_exchange>("plane", nplanes);

			sman.set<services::message_passing>("plane-setup");
			sman.set<services::mempool>("plane");
			if (aggregate)
				sman.set<mpi::halo_exchange>("plane");
		} else {
			service_manager<cdr2::mpi::stypes> & master = kmans[0]->services();
			auto mp_service = master.get_ptr<services::message_passing>();
			auto mp_service_solve = master.get_ptr<services::message_passing>("plane");
			auto mempool_service = master.get_ptr<services::mempool>();
			auto halo_service = master.get_ptr<services::halo_exchange<cdr2::mpi::stypes>>("plane");
			auto *mpi_keys = static_cast<cdr3::mpi::plane_setup_mpi*>(mp_service.get())->get_keys();
			auto *addrs = static_cast<cdr3::mpi::plane_mempool*>(mempool_service.get())->get_addrs();
			auto barrier = static_cast<cdr3::mpi::plane_exchange*>(halo_service.get())->get_barrier();
			auto barrier_mpi = static_cast<cdr3::mpi::plane_mpi*>(mp_service_solve.get())->get_barrier();
			sman.add<services::message_passing, cdr3::mpi::plane_setup_mpi>("plane-setup", mpi_keys);
			sman.add<services::message_passing, cdr3::mpi::plane_mpi>("plane", nplanes, barrier_mpi);
			sman.add<services::mempool, cdr3::mpi::plane_mempool>("plane", i, addrs);
			sman.add<mpi::halo_exchange, cdr3::mpi::plane_exchange>("plane", nplanes, barrier);

			sman.set<services::message_passing>("plane-setup");
			sman.set<services::mempool>("plane");
			if (aggregate)
				sman.set<mpi::halo_exchange>("plane");
		}
	}

	for (auto i : range<std::size_t>(nplanes)) {
		auto & sman = kmans[i]->services();
		auto & halo = sman.get<mpi::halo_exchange>();
		auto & mpool = sman.get<services::mempool>();

		halo.setup({{grid}});

		std::size_t nbytes = grid->nlocal(0) * grid->nlocal(1) * sizeof(real_t);
		real_t *xaddr = (real_t*) mpool.addr(services::mempool::sol, nbytes);
		real_t *baddr = (real_t*) mpool.addr(services::mempool::rhs, nbytes);

		xs.emplace_back(xaddr, grid);
		bs.emplace_back(baddr, grid);

		xs[i].set(0.0);
		bs[i].set(i+1);
	}

	auto so = mpi::gallery::poisson(grid);
	mpi::grid_func res(grid);
	relax_stencil sor(grid->nlocal(0) - 2, grid->nlocal(1) - 2);

	for (auto i : range<std::size_t>(nplanes)) {
		kmans[i]->setup<mpi::line_relax<relax_dir::x>>(so, sor);
	}


	if (aggregate) {
		for (auto i : range<std::size_t>(nplanes)) {
			auto & sman = kmans[i]->services();
			sman.set<services::message_passing>("plane");
		}
	}

	ABT_thread *threads;
	if (aggregate) {
		threads = new ABT_thread[nplanes];
		ult_params *args = new ult_params[nplanes];

		for (auto i : range<std::size_t>(nplanes)) {
			args[i].A = &so;
			args[i].x = &xs[i];
			args[i].b = &bs[i];
			args[i].sor = &sor;
			args[i].res = &res;
			args[i].kman = kmans[i];
			ABT_thread_create(pool, compute, (void*) &args[i],
			                  ABT_THREAD_ATTR_NULL, &threads[i]);
		}

		for (auto i : range<std::size_t>(nplanes)) {
			ABT_thread_join(threads[i]);
		}

		for (auto i : range<std::size_t>(nplanes)) {
			ABT_thread_free(&threads[i]);
		}

		delete[] threads;
		delete[] args;
	} else {
		for (auto i : range<std::size_t>(nplanes)) {
			kmans[i]->run<mpi::line_relax<relax_dir::x>>(so, xs[i], bs[i], sor, res, cycle::Dir::DOWN);
		}
	}

	// int rank;
	// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// if (rank == 0) {
	// 	std::cout << xs[1] << std::endl;
	// }

	ABT_finalize();
	MPI_Finalize();
	return 0;
}
