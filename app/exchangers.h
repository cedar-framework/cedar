#ifndef TEST_EXCHANGERS_H
#define TEST_EXCHANGERS_H

#include <cedar/2d/mpi/msg_exchanger.h>
#include <cedar/2d/mpi/tausch_exchanger.h>

namespace cedar { namespace cdr2 { namespace mpi {
class noop_exchange : public msg_exchanger
{
public:
	void setup(std::vector<topo_ptr> topos) override {}
	void run(stencil_op<five_pt> & so) override {}
	void run(stencil_op<nine_pt> & so) override {}
	void run(grid_func & gf) override {
		#pragma omp barrier
		#pragma omp barrier
	}
	void exchange_func(int k, real_t *gf) override {
		#pragma omp barrier
		#pragma omp barrier
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
		#pragma omp barrier
		parent::run(gf);
		#pragma omp barrier
	}


	void exchange_func(int k, real_t *gf) override {
		#pragma omp barrier
		parent::exchange_func(k, gf);
		#pragma omp barrier
	}

private:
	int nplanes;
};

}}}

#endif
