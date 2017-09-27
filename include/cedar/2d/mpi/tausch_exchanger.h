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
	tausch_exchanger(const kernel_params & params, grid_topo & topo);
	virtual void exchange_func(int k, real_t * gf);
	virtual void exchange_sten(int k, real_t * so);
	template<class sten>
		void exchange(mpi::stencil_op<sten> & so);
	void exchange(mpi::grid_func & f);

protected:
	std::unique_ptr<Tausch<double>> tausch;
	std::array<bool, halo_dir::count> active;
};

template<class sten>
	void tausch_exchanger::exchange(mpi::stencil_op<sten> & so)
{

}

}}}


#endif
