#ifndef CEDAR_TEST_LINE_AGG_H
#define CEDAR_TEST_LINE_AGG_H

#include <abt.h>

#include <cedar/2d/mpi/kernel_manager.h>
#include <cedar/2d/mpi/gallery.h>
#include <cedar/2d/util/topo.h>

#include <cedar/3d/mpi/plane_mempool.h>
#include <cedar/3d/mpi/plane_mpi.h>
#include <cedar/3d/mpi/plane_exchange.h>

using namespace cedar;
using namespace cedar::cdr2;

template<class sten>
mpi::stencil_op<sten> create_op(cedar::topo_ptr grid);

template<class sten2>
struct ult_params
{
	mpi::stencil_op<sten2> *A;
	mpi::grid_func *x;
	mpi::grid_func *b;
	relax_stencil *sor;
	mpi::grid_func *res;
	mpi::kman_ptr kman;
	int nsweeps;
};

template<cedar::relax_dir rdir, class sten2>
void compute(void *args)
{
	ult_params<sten2> *params = (ult_params<sten2>*) args;

	auto kman = params->kman;
	int nsweeps = params->nsweeps;
	auto & so = *params->A;
	auto & x = *params->x;
	auto & b = *params->b;
	auto & sor = *params->sor;
	auto & res = *params->res;

	for (auto sweep : range(nsweeps)) {
		kman->template run<mpi::line_relax<rdir>>(so, x, b, sor, res, cycle::Dir::DOWN);
		kman->template run<mpi::line_relax<rdir>>(so, x, b, sor, res, cycle::Dir::UP);
	}
}


template<cedar::relax_dir rdir, class sten2>
std::vector<real_t> line_agg(bool aggregate, int nx, int ny, int nplanes)
{
	auto grid = util::create_topo_global(MPI_COMM_WORLD, nx, ny);
	config conf("test-planeagg-lines.json");
	log::init(conf);
	int nsweeps = 1;

	ABT_pool pool;
	ABT_xstream xstream;
	ABT_xstream_self(&xstream);
	ABT_xstream_get_main_pools(xstream, 1, &pool);

	std::vector<ABT_thread> threads;
	threads.reserve(nplanes);

	std::vector<mpi::kman_ptr> kmans;
	std::vector<mpi::grid_func> bs;
	std::vector<mpi::grid_func> xs;
	std::vector<mpi::grid_func> res;

	for (auto i : range<std::size_t>(nplanes)) {
		kmans.push_back(mpi::build_kernel_manager(conf));
		auto & sman = kmans.back()->services();
		if (i == 0) {
			sman.add<services::message_passing, cdr3::mpi::plane_setup_mpi>("plane-setup");
			sman.add<services::message_passing, cdr3::mpi::plane_mpi>("plane", nplanes, &threads);
			sman.add<services::mempool, cdr3::mpi::plane_mempool>("plane", nplanes);
			sman.add<mpi::halo_exchange, cdr3::mpi::plane_exchange>("plane", nplanes, &threads);

			sman.set<services::message_passing>("plane-setup");
			sman.set<services::mempool>("plane");
			if (aggregate)
				sman.set<mpi::halo_exchange>("plane");
		} else {
			service_manager<cdr2::mpi::stypes> & master = kmans[0]->services();
			auto mp_service = master.get_ptr<services::message_passing>();
			auto mempool_service = master.get_ptr<services::mempool>();
			auto *mpi_keys = static_cast<cdr3::mpi::plane_setup_mpi*>(mp_service.get())->get_keys();
			auto *addrs = static_cast<cdr3::mpi::plane_mempool*>(mempool_service.get())->get_addrs();

			sman.add<services::message_passing, cdr3::mpi::plane_setup_mpi>("plane-setup", mpi_keys);
			sman.add<services::message_passing, cdr3::mpi::plane_mpi>("plane", nplanes, i, &threads);
			sman.add<services::mempool, cdr3::mpi::plane_mempool>("plane", i, addrs);
			sman.add<mpi::halo_exchange, cdr3::mpi::plane_exchange>("plane", nplanes, i, &threads);

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
		real_t *raddr = (real_t*) mpool.addr(services::mempool::res, nbytes);

		xs.emplace_back(xaddr, grid);
		bs.emplace_back(baddr, grid);
		res.emplace_back(raddr, grid);

		xs[i].set(0.0);
		bs[i].set(i+1);
	}

	auto so = create_op<sten2>(grid);
	relax_stencil sor(grid->nlocal(0) - 2, grid->nlocal(1) - 2);

	for (auto i : range<std::size_t>(nplanes)) {
		kmans[i]->setup<mpi::line_relax<rdir>>(so, sor);
	}


	if (aggregate) {
		for (auto i : range<std::size_t>(nplanes)) {
			auto & sman = kmans[i]->services();
			sman.set<services::message_passing>("plane");
		}
	}

	if (aggregate) {
		ult_params<sten2> *args = new ult_params<sten2>[nplanes];

		for (auto i : range<std::size_t>(nplanes)) {
			args[i].A = &so;
			args[i].x = &xs[i];
			args[i].b = &bs[i];
			args[i].sor = &sor;
			args[i].res = &res[i];
			args[i].kman = kmans[i];
			args[i].nsweeps = nsweeps;
		}

		for (auto i : range<std::size_t>(nplanes)) {
			ABT_thread_create(pool, compute<rdir, sten2>, (void*) &args[i],
			                  ABT_THREAD_ATTR_NULL, &threads[i]);
		}

		for (auto i : range<std::size_t>(nplanes)) {
			ABT_thread_join(threads[i]);
		}

		for (auto i : range<std::size_t>(nplanes)) {
			ABT_thread_free(&threads[i]);
		}

		delete[] args;
	} else {
		for (auto sweep : range(nsweeps)) {
			for (auto i : range<std::size_t>(nplanes)) {
				kmans[i]->run<mpi::line_relax<rdir>>(so, xs[i], bs[i], sor, res[i], cycle::Dir::DOWN);
				kmans[i]->run<mpi::line_relax<rdir>>(so, xs[i], bs[i], sor, res[i], cycle::Dir::UP);
			}
		}
	}

	// need to copy since mempool will free memory when it loses scope
	std::vector<real_t> ret;
	for (auto & xcur : xs) {
		for (auto j : xcur.range(1)) {
			for (auto i : xcur.range(0)) {
				ret.push_back(xcur(i,j));
			}
		}
	}

	return ret;
}

#endif
