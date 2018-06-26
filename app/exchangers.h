#ifndef TEST_EXCHANGERS_H
#define TEST_EXCHANGERS_H

#include <abt.h>
#include <cedar/2d/mpi/msg_exchanger.h>
#include <cedar/2d/mpi/tausch_exchanger.h>

ABT_barrier comm_barrier;

namespace cedar { namespace cdr2 { namespace mpi {
class noop_exchange : public msg_exchanger
{
public:
	void setup(std::vector<topo_ptr> topos) override {}
	void run(stencil_op<five_pt> & so) override {}
	void run(stencil_op<nine_pt> & so) override {}
	void run(grid_func & gf) override {
		ABT_barrier_wait(comm_barrier);
		ABT_barrier_wait(comm_barrier);
	}
	void exchange_func(int k, real_t *gf) override {
		ABT_barrier_wait(comm_barrier);
		ABT_barrier_wait(comm_barrier);
	}
	void exchange_sten(int k, real_t * so) override {}
};


class master_exchange : public tausch_exchanger
{
public:
	using parent = tausch_exchanger;
	master_exchange(int nplanes): nplanes(nplanes) {}
	void setup(std::vector<topo_ptr> topos) override
	{
		ABT_barrier_create((std::size_t) nplanes, &comm_barrier);
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


	void run(grid_func & gf) override {
		ABT_barrier_wait(comm_barrier);
		parent::run(gf);
		ABT_barrier_wait(comm_barrier);
	}


	void exchange_func(int k, real_t *gf) override {
		ABT_barrier_wait(comm_barrier);
		parent::exchange_func(k, gf);
		ABT_barrier_wait(comm_barrier);
	}

private:
	int nplanes;
};

}}}

#endif
