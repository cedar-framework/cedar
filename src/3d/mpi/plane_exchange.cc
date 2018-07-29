#include <cedar/3d/mpi/plane_exchange.h>

using namespace cedar;
using namespace cedar::cdr3::mpi;


plane_exchange::plane_exchange(int nplanes) :
	nplanes(nplanes), ismaster(true)
{
	if (ABT_initialized() == ABT_ERR_UNINITIALIZED)
			ABT_init(0, NULL);

	ABT_barrier_create((std::size_t) nplanes, &barrier);
	halo3 = std::make_unique<tausch_exchanger>();
}


plane_exchange::plane_exchange(int nplanes, ABT_barrier barrier) :
	nplanes(nplanes), ismaster(false), barrier(barrier) {
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
	halo3->setup(topos);
	halo3->activate_send(cdr3::mpi::tausch_exchanger::halo_dir::top, false);
	halo3->activate_recv(cdr3::mpi::tausch_exchanger::halo_dir::bottom, false);

	for (auto tpr : topos) {
		tpr->nlocal(2) = 0;
	}
}


void plane_exchange::run(grid_func & gf)
{
	ABT_barrier_wait(barrier);
	if (ismaster) {
		auto topo = gf.grid_ptr();
		cdr3::mpi::grid_func gf3(gf.data(), topo);
		halo3->run(gf3);
	}
	ABT_barrier_wait(barrier);
}


void plane_exchange::exchange_func(int k, real_t *gf)
{
	ABT_barrier_wait(barrier);
	if (ismaster)
		halo3->exchange_func(k, gf);
	ABT_barrier_wait(barrier);
}
