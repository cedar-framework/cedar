#ifndef CEDAR_2D_TAUSCH_EXCHANGER_H
#define CEDAR_2D_TAUSCH_EXCHANGER_H

#include <tausch/tausch.h>

#include <cedar/kernels/halo_exchange.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/mpi/types.h>
#include <cedar/mpi/grid_topo.h>

namespace cedar { namespace cdr2 { namespace mpi {


struct line_pkg
{
	line_pkg(grid_topo & topo);
	std::vector<real_t> linebuf;
	std::array<MPI_Comm, 2> linecomm;
	std::array<array<len_t, 2>, 2> datadist;
};

class tausch_exchanger : public kernels::halo_exchange<stypes>
{
public:
	enum halo_dir { left=0, right=1, down=2, up=3, count };
	void setup(std::vector<topo_ptr> topos) override;
	void run(stencil_op<five_pt> & so) override { this->run_impl(so); }
	void run(stencil_op<nine_pt> & so) override { this->run_impl(so); }
	void run(grid_func & gf) override;
	void exchange_func(int k, real_t * gf) override;
	void exchange_sten(int k, real_t * so) override;
	aarray<int, len_t, 2> & leveldims(int k) override {
		return dims[k];
	}
	len_t * datadist(int k, int grid) override {
		return line_data->datadist[k].data();
	}
	MPI_Comm linecomm(int k) override {
		return line_data->linecomm[k];
	}
	std::vector<real_t> & linebuf() override {
		return line_data->linebuf;
	}

	template<class sten>
	void run_impl(stencil_op<sten> & so)
	{
		auto lvl = so.grid().level();
		lvl = nlevels - lvl - 1;

		for (int dir = 0; dir < halo_dir::count; dir++) {
			if (recv_active[index(lvl, dir)])
				tausch_so->postReceive(TAUSCH_CwC, index(lvl, dir), index(lvl, dir));
		}

		for (int dir = 0; dir < halo_dir::count; dir++) {
			if (send_active[index(lvl, dir)]) {
				for (int sdir = 0; sdir < stencil_ndirs<sten>::value; sdir++)
					tausch_so->packSendBuffer(TAUSCH_CwC, index(lvl,dir), sdir, so.data() + so.index(0,0,sdir));
				tausch_so->send(TAUSCH_CwC, index(lvl,dir), index(lvl,dir));

			}
			if (recv_active[index(lvl, dir)]) {
				tausch_so->recv(TAUSCH_CwC, index(lvl,dir));
				for (int sdir = 0; sdir < stencil_ndirs<sten>::value; sdir++)
					tausch_so->unpackRecvBuffer(TAUSCH_CwC, index(lvl,dir), sdir, so.data() + so.index(0,0,sdir));
			}
		}
	}
	std::unique_ptr<line_pkg> line_data;

protected:
	std::unique_ptr<Tausch<real_t>> tausch;
	std::unique_ptr<Tausch<real_t>> tausch_so;
	std::vector<bool> send_active;
	std::vector<bool> recv_active;

	void set_level_spec(int lvl, int rank,
	                    grid_topo & topo,
	                    std::vector<TauschHaloSpec> & remote_spec,
	                    std::vector<TauschHaloSpec> & local_spec);
	void set_level_spec_so(int lvl, int rank,
	                       grid_topo & topo,
	                       std::vector<TauschHaloSpec> & remote_spec,
	                       std::vector<TauschHaloSpec> & local_spec);
	void init_gfunc(std::vector<topo_ptr> & topos);
	void init_so(std::vector<topo_ptr> & topos);
	void init_dims(grid_topo & topo);
	void init_datadist();
	std::size_t index(int lvl, int dir) { return lvl*halo_dir::count + dir; }
	std::size_t nlevels;
	std::array<aarray<int, len_t, 2>, 2> dims;
	std::array<std::vector<len_t>, 2> dimfine;
	std::array<int, 2> coord;
};

}}}


#endif
