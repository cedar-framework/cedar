#ifndef TEST_EXCHANGERS_H
#define TEST_EXCHANGERS_H

#include <abt.h>
#include <cedar/2d/mpi/msg_exchanger.h>
#include <cedar/2d/mpi/tausch_exchanger.h>
#include <cedar/3d/mpi/tausch_exchanger.h>

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


class master_exchange : public services::halo_exchange<stypes>
{
public:
	master_exchange(int nplanes, bool reimpl): nplanes(nplanes), reimpl(reimpl) {
		halo3 = std::make_unique<cdr3::mpi::tausch_exchanger>();
	}
	void setup(std::vector<topo_ptr> topos) override
	{
		std::vector<topo_ptr> topos_all;

		for (auto tpr : topos) {
			tpr->nlocal(2) = nplanes;
			tpr->coord(2) = 0;
			tpr->nproc(2) = 1;
		}

		halo3->add_params(this->params);
		halo3->setup(topos);

		for (auto tpr : topos) {
			tpr->nlocal(2) = 1;
		}

		halo3->activate_send(cdr3::mpi::tausch_exchanger::halo_dir::top, false);
		halo3->activate_recv(cdr3::mpi::tausch_exchanger::halo_dir::bottom, false);
	}


	void run(grid_func & gf) override {
		auto topo = gf.grid_ptr();
		cdr3::mpi::grid_func gf3(gf.data(), topo);
		halo3->run(gf3);
		if (not reimpl)
			ABT_thread_yield_to(threads[0]);
	}


	void exchange_func(int k, real_t *gf) override {
		halo3->exchange_func(k, gf);
		if (not reimpl)
			ABT_thread_yield_to(threads[0]);
	}

	aarray<int, len_t, 2> & leveldims(int k) override { return halo3->leveldims(k); }
	len_t * datadist(int k, int grid) override { return halo3->datadist(k, grid); }
	std::vector<real_t> & linebuf() override { return halo3->linebuf(); }

	void run(mpi::stencil_op<nine_pt> & so) override {}
	void run(mpi::stencil_op<five_pt> & so) override {}
	void exchange_sten(int k, real_t *so) override {}

private:
	std::unique_ptr<cdr3::mpi::tausch_exchanger> halo3;
	int nplanes;
	bool reimpl;
};

}}}

#endif
