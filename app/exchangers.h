#ifndef TEST_EXCHANGERS_H
#define TEST_EXCHANGERS_H

#include <abt.h>
#include <cedar/2d/mpi/msg_exchanger.h>
#include <cedar/2d/mpi/tausch_exchanger.h>

ABT_barrier comm_barrier;
ABT_thread *threads;

namespace cedar { namespace cdr2 { namespace mpi {
class noop_exchange : public msg_exchanger
{
public:
	void setup(std::vector<topo_ptr> topos) override {}
	void run(stencil_op<five_pt> & so) override {}
	void run(stencil_op<nine_pt> & so) override {}
	void run(grid_func & gf) override {
		ABT_thread_id tid;
		ABT_thread_self_id(&tid);
		ABT_thread_yield_to(threads[tid+1]);
	}
	void exchange_func(int k, real_t *gf) override {
		ABT_thread_id tid;
		ABT_thread_self_id(&tid);
		ABT_thread_yield_to(threads[tid+1]);
	}
	void exchange_sten(int k, real_t * so) override {}
};


class master_exchange : public tausch_exchanger
{
public:
	using parent = tausch_exchanger;
	master_exchange(int nplanes, bool reimpl): nplanes(nplanes), reimpl(reimpl) {}
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
		parent::run(gf);
		if (not reimpl)
			ABT_thread_yield_to(threads[0]);
	}


	void exchange_func(int k, real_t *gf) override {
		parent::exchange_func(k, gf);
		if (not reimpl)
			ABT_thread_yield_to(threads[0]);
	}

private:
	int nplanes;
	bool reimpl;
};

}}}

#endif
