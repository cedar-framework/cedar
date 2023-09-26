#ifndef CEDAR_2D_MPI_RESTRICT_H
#define CEDAR_2D_MPI_RESTRICT_H

#include <cedar/kernels/restrict.h>
#include <cedar/2d/mpi/types.h>
#include <cedar/device.h>

namespace cedar { namespace cdr2 { namespace mpi {

template <typename device=cedar::cpu>
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

}}}

#endif
