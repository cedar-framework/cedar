#ifndef CEDAR_KERNELS_SETUP_NOG_H
#define CEDAR_KERNELS_SETUP_NOG_H

#include <cedar/kernel.h>
#include <cedar/mpi/grid_topo.h>

namespace cedar { namespace kernels {

template<class solver_types>
class setup_nog : public kernel<solver_types>
{
public:
	const std::string name = "setup nog";

	virtual void run(grid_topo & topo,
	                 len_t min_coarse, int *nog) = 0;
};

}}
#endif

