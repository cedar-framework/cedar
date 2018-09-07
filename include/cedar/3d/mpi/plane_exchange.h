#ifndef CEDAR_3D_MPI_PLANE_EXCHANGE_H
#define CEDAR_3D_MPI_PLANE_EXCHANGE_H

#include <abt.h>

#include <cedar/2d/mpi/types.h>
#include <cedar/3d/mpi/tausch_exchanger.h>

namespace cedar { namespace cdr3 { namespace mpi {

class plane_exchange : public services::halo_exchange<cdr2::mpi::stypes>
{
public:
	plane_exchange(int nplanes, std::vector<ABT_thread> *threads);
	plane_exchange(int nplanes, int worker_id, std::vector<ABT_thread> *threads);

	void setup(std::vector<topo_ptr> topos) override;
	void run(stencil_op<cdr2::five_pt> & so) override
		{ log::error << "plane exchange of stencil_op not supported" << std::endl; }
	void run(stencil_op<cdr2::nine_pt> & so) override
		{ log::error << "plane exchange of stencil_op not supported" << std::endl; }
	void run(grid_func & gf, unsigned short dmask=7) override;
	void exchange_func(int k, real_t *gf) override;
	void exchange_sten(int k, real_t *gf) override
		{ log::error << "plane exchange of stencil_op not supported" << std::endl; }

	aarray<int, len_t, 2> & leveldims(int k) override { return halo3->leveldims(k); }
	len_t * datadist(int k, int grid) override { return halo3->datadist(k, grid); }
	std::vector<real_t> & linebuf() override { return halo3->linebuf(); }

protected:
	std::unique_ptr<tausch_exchanger> halo3;
	int nplanes;
	bool ismaster;
	std::vector<ABT_thread> *threads;
	int wid;
};

}}}

#endif
