#ifndef CEDAR_2D_TAUSCH_EXCHANGER_H
#define CEDAR_2D_TAUSCH_EXCHANGER_H

#include <tausch/tausch.h>

#include <cedar/kernel_params.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/halo_exchanger_base.h>

namespace cedar { namespace cdr2 { namespace mpi {

class tausch_exchanger : public halo_exchanger_base
{
	enum halo_dir { left=0, right=1, down=2, up=3, count };
public:
	tausch_exchanger(const kernel_params & params,
	                 std::vector<topo_ptr> topo);
	virtual void exchange_func(int k, real_t * gf);
	virtual void exchange_sten(int k, real_t * so);
	template<class sten>
		void exchange(mpi::stencil_op<sten> & so);
	void exchange(mpi::grid_func & f);
	virtual aarray<int, len_t, 2> & leveldims(int k) {
		return dims[k];
	}

protected:
	std::unique_ptr<Tausch<real_t>> tausch;
	std::unique_ptr<Tausch<real_t>> tausch_so;
	std::vector<bool> send_active;
	std::vector<bool> recv_active;

private:
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
	std::size_t index(int lvl, int dir) { return lvl*halo_dir::count + dir; }
	std::size_t nlevels;
	/* aarray<int, len_t, 2> dimx; */
	/* aarray<int, len_t, 2> dimy; */
	std::array<aarray<int, len_t, 2>, 2> dims;
	std::array<std::vector<len_t>, 2> dimfine;
};

template<class sten>
	void tausch_exchanger::exchange(mpi::stencil_op<sten> & so)
{
	auto lvl = so.grid().level();

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

}}}


#endif
