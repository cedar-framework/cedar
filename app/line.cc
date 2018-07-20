#include <cedar/2d/mpi/kernel_manager.h>
#include <cedar/2d/mpi/gallery.h>
#include <cedar/2d/util/topo.h>

#include <cedar/3d/mpi/plane_mempool.h>
#include <cedar/3d/mpi/plane_mpi.h>
#include <cedar/3d/mpi/plane_exchange.h>

int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr2;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);

	config conf("config.json");
	log::init(conf);

	auto nplanes = conf.get<int>("nplanes");

	auto grid = util::create_topo(conf);

	std::vector<mpi::kman_ptr> kmans;
	std::vector<mpi::grid_func> bs;
	std::vector<mpi::grid_func> xs;

	for (auto i : range<std::size_t>(nplanes)) {
		kmans.push_back(mpi::build_kernel_manager(conf));
		auto & sman = kmans.back()->services();
		if (i == 0) {
			sman.add<services::message_passing, cdr3::mpi::plane_setup_mpi>("plane-setup");
			sman.add<services::mempool, cdr3::mpi::plane_mempool>("plane", nplanes);
			sman.add<services::halo_exchange<cdr2::mpi::stypes>, cdr3::mpi::plane_exchange>("plane", nplanes);

			sman.set<services::message_passing>("plane-setup");
			sman.set<services::mempool>("plane");
			sman.set<services::halo_exchange<cdr2::mpi::stypes>>("plane");
		} else {
			service_manager<cdr2::mpi::stypes> & master = kmans[0]->services();
			auto mp_service = master.get_ptr<services::message_passing>();
			auto mempool_service = master.get_ptr<services::mempool>();
			auto halo_service = master.get_ptr<services::halo_exchange<cdr2::mpi::stypes>>("plane");
			auto *mpi_keys = static_cast<cdr3::mpi::plane_setup_mpi*>(mp_service.get())->get_keys();
			auto *addrs = static_cast<cdr3::mpi::plane_mempool*>(mempool_service.get())->get_addrs();
			auto barrier = static_cast<cdr3::mpi::plane_exchange*>(halo_service.get())->get_barrier();
			sman.add<services::message_passing, cdr3::mpi::plane_setup_mpi>("plane-setup", mpi_keys);
			sman.add<services::mempool, cdr3::mpi::plane_mempool>("plane", i, addrs);
			sman.add<services::halo_exchange<cdr2::mpi::stypes>, cdr3::mpi::plane_exchange>("plane", nplanes, barrier);

			sman.set<services::message_passing>("plane-setup");
			sman.set<services::mempool>("plane");
			sman.set<services::halo_exchange<cdr2::mpi::stypes>>("plane");
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
		bs[i].set(i);
	}

	auto so = mpi::gallery::poisson(grid);
	mpi::grid_func res(grid);
	relax_stencil sor(grid->nlocal(0) - 2, grid->nlocal(1) - 2);

	for (auto i : range<std::size_t>(nplanes)) {
		kmans[i]->setup<mpi::line_relax<relax_dir::x>>(so, sor);
	}

	for (auto i : range<std::size_t>(nplanes)) {
		kmans[i]->run<mpi::line_relax<relax_dir::x>>(so, xs[i], bs[i], sor, res, cycle::Dir::DOWN);
	}

	MPI_Finalize();
	return 0;
}
