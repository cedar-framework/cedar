#include <cedar/3d/mpi/plane_exchange.h>

using namespace cedar;
using namespace cedar::cdr3::mpi;


void plane_exchange::setup(std::vector<topo_ptr> topos)
{
	if (ismaster) {
		if (ABT_initialized() == ABT_ERR_UNINITIALIZED)
			ABT_init(0, NULL);

		ABT_barrier_create((std::size_t) nplanes, &barrier);
		std::vector<topo_ptr> topos_all;

		for (auto tpr : topos) {
			tpr->nlocal(1) *= nplanes;
		}

		parent::setup(topos);

		for (auto tpr : topos) {
			tpr->nlocal(1) /= nplanes;
		}

		for (std::size_t lvl = 0; lvl < nlevels; lvl++) {
			send_active[index(lvl, halo_dir::down)] = false;
			recv_active[index(lvl, halo_dir::down)] = false;
			send_active[index(lvl, halo_dir::up)] = false;
			recv_active[index(lvl, halo_dir::up)] = false;
		}
	}
}


void plane_exchange::run(grid_func & gf)
{
	ABT_barrier_wait(barrier);
	if (ismaster)
		parent::run(gf);
	ABT_barrier_wait(barrier);
}


void plane_exchange::exchange_func(int k, real_t *gf)
{
	ABT_barrier_wait(barrier);
	if (ismaster)
		parent::exchange_func(k, gf);
	ABT_barrier_wait(barrier);
}
