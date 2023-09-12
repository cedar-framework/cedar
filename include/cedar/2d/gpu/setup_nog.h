#ifndef CEDAR_2D_MPI_SETUP_NOG_H
#define CEDAR_2D_MPI_SETUP_NOG_H

#include <cedar/2d/gpu/types.h>
#include <cedar/kernels/setup_nog.h>

namespace cedar { namespace cdr2 { namespace gpu { namespace mpi {

class setup_nog_f90 : public kernels::setup_nog<stypes>
{
	void run(grid_topo & topo, len_t min_coarse, int *nog) override;
};

}}}}

#endif
