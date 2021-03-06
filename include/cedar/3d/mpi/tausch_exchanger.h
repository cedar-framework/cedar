#ifndef CEDAR_3D_TAUSCH_EXCHANGER_H
#define CEDAR_3D_TAUSCH_EXCHANGER_H

#include <tausch.h>

#include <cedar/services/halo_exchange.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/grid_func.h>
#include <cedar/3d/mpi/types.h>
#include <cedar/mpi/grid_topo.h>

namespace cedar { namespace cdr3 { namespace mpi {

struct line_pkg
{
	line_pkg(grid_topo & topo);
	std::vector<real_t> linebuf;
	std::array<array<len_t, 2>, 3> datadist;
};


class tausch_exchanger : public services::halo_exchange<stypes>
{
public:
	using parent = services::halo_exchange<stypes>;
	enum halo_dir { west=0, east=1, south=2, north=3, bottom=4, top=5, count };
	tausch_exchanger() : nghost{{2,2,2}} {}
	void setup_worker(std::vector<topo_ptr> topos);
	void setup(std::vector<topo_ptr> topos) override;
	void run(stencil_op<seven_pt> & so) override { this->run_impl(so); }
	void run(stencil_op<xxvii_pt> & so) override { this->run_impl(so); }
	void run(grid_func & gf, unsigned short dmask=7) override;
	void exchange_func(int k, real_t * gf) override;
	void exchange_sten(int k, real_t * so) override;
	aarray<int, len_t, 2> & leveldims(int k) override {
		return dims[k];
	}
	len_t * datadist(int k, int grid) override {
		return line_data->datadist[k].data();
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
			if (send_active[index(lvl, dir)]) {
				for (int sdir = 0; sdir < stencil_ndirs<sten>::value; sdir++)
					tausch_so->packSendBuffer(index(lvl,dir), sdir, so.data() + so.index(0,0,sdir));
				tausch_so->send(index(lvl,dir), index(lvl,dir));

			}
			if (recv_active[index(lvl, dir)]) {
				tausch_so->recv(index(lvl,dir), index(lvl,dir));
				for (int sdir = 0; sdir < stencil_ndirs<sten>::value; sdir++)
					tausch_so->unpackRecvBuffer(index(lvl,dir), sdir, so.data() + so.index(0,0,sdir));
			}
		}
	}
	std::unique_ptr<line_pkg> line_data;

	void activate_send(halo_dir dir, bool active);
	void activate_recv(halo_dir dir, bool active);

	int & nhalo(int dim) { return nghost[dim]; }

protected:
	std::unique_ptr<Tausch> tausch;
	std::unique_ptr<Tausch> tausch_so;
	std::vector<bool> send_active;
	std::vector<bool> recv_active;

	void set_level_spec(int lvl, int rank,
	                    grid_topo & topo,
	                    std::vector<std::vector<std::array<int, 4> > > & remote_spec,
	                    std::vector<std::vector<std::array<int, 4> > > & local_spec,
	                    std::vector<int> & remote_remoteMpiRank,
	                    std::vector<int> & local_remoteMpiRank);
	void set_level_spec_so(int lvl, int rank,
	                       grid_topo & topo,
	                       std::vector<std::vector<std::array<int, 4> > > & remote_spec,
	                       std::vector<std::vector<std::array<int, 4> > > & local_spec,
	                       std::vector<int> & remote_remoteMpiRank,
	                       std::vector<int> & local_remoteMpiRank);
	void init_gfunc(std::vector<topo_ptr> & topos);
	void init_so(std::vector<topo_ptr> & topos);
	void init_dims(grid_topo & topo);
	void init_datadist();
	std::size_t index(int lvl, int dir) { return lvl*halo_dir::count + dir; }
	std::size_t nlevels;
	std::array<aarray<int, len_t, 2>, 3> dims;
	std::array<std::vector<len_t>, 3> dimfine;
	std::array<int, 3> coord;
	std::array<int, 3> nghost;
};


}}}
#endif
