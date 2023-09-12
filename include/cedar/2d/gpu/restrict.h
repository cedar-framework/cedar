#ifndef CEDAR_2D_GPU_RESTRICT_H
#define CEDAR_2D_GPU_RESTRICT_H

#include <cedar/kernels/restrict.h>
#include <cedar/2d/gpu/types.h>

namespace cedar { namespace cdr2 { namespace gpu { namespace mpi {

class restrict_f90 : public kernels::restriction<stypes>
{
	void run(const restrict_op & R,
	         const grid_func & x,
	         grid_func & y) override;

	void update_periodic(mpi::grid_func & q,
	                     const grid_topo & topof,
	                     const grid_topo & topoc,
	                     std::array<bool, 3> periodic);
};

}}}}

#endif
