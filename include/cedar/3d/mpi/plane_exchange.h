#ifndef CEDAR_3D_MPI_PLANE_EXCHANGE_H
#define CEDAR_3D_MPI_PLANE_EXCHANGE_H

#include <abt.h>

#include <cedar/2d/mpi/tausch_exchanger.h>

namespace cedar { namespace cdr3 { namespace mpi {

class plane_exchange : public cdr2::mpi::tausch_exchanger
{
public:
	static ABT_barrier barrier;
	using parent = cdr2::mpi::tausch_exchanger;

	plane_exchange(int nplanes, bool ismaster) :
		nplanes(nplanes), ismaster(ismaster) {}

	void setup(std::vector<topo_ptr> topos) override;
	void run(stencil_op<cdr2::five_pt> & so) override
		{ log::error << "plane exchange of stencil_op not supported" << std::endl; }
	void run(stencil_op<cdr2::nine_pt> & so) override
		{ log::error << "plane exchange of stencil_op not supported" << std::endl; }
	void run(grid_func & gf) override;
	void exchange_func(int k, real_t *gf) override;
	void exchange_sten(int k, real_t *gf) override
		{ log::error << "plane exchange of stencil_op not supported" << std::endl; }

protected:
	int nplanes;
	bool ismaster;
};

}}}

#endif
