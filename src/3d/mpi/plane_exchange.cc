#include <cedar/3d/mpi/plane_exchange.h>

using namespace cedar;
using namespace cedar::cdr3::mpi;


plane_exchange::plane_exchange(int nplanes, std::vector<ABT_thread> *threads) :
	nplanes(nplanes), ismaster(true), threads(threads), wid(0)
{
	if (ABT_initialized() == ABT_ERR_UNINITIALIZED)
			ABT_init(0, NULL);

	halo3 = std::make_unique<tausch_exchanger>();
}


plane_exchange::plane_exchange(int nplanes, int worker_id, std::vector<ABT_thread> *threads) :
	nplanes(nplanes), ismaster(false), threads(threads), wid(worker_id) {
	halo3 = std::make_unique<tausch_exchanger>();
}


void plane_exchange::setup(std::vector<topo_ptr> topos)
{
	for (auto tpr : topos) {
		tpr->nlocal(2) = nplanes;
		tpr->coord(2) = 0;
		tpr->nproc(2) = 1;
	}

	halo3->nhalo(2) = 0;
	halo3->add_params(this->params);
	if (ismaster) {
		halo3->setup(topos);
		halo3->activate_send(cdr3::mpi::tausch_exchanger::halo_dir::top, false);
		halo3->activate_recv(cdr3::mpi::tausch_exchanger::halo_dir::bottom, false);
	} else {
		halo3->setup_worker(topos);
	}

	for (auto tpr : topos) {
		tpr->nlocal(2) = 0;
	}
}


void plane_exchange::run(grid_func & gf, unsigned short dmask)
{
	ABT_thread_yield_to((*threads)[(wid+1) % nplanes]);
	if (ismaster) {
		auto topo = gf.grid_ptr();
		cdr3::mpi::grid_func gf3(gf.data(), topo);
		halo3->run(gf3);
	}
}


void plane_exchange::exchange_func(int k, real_t *gf)
{
	ABT_thread_yield_to((*threads)[(wid+1) % nplanes]);
	if (ismaster)
		halo3->exchange_func(k, gf);
}
