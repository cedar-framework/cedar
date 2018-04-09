#ifndef CEDAR_3D_MPI_RESTRICT_H
#define CEDAR_3D_MPI_RESTRICT_H

#include <cedar/3d/mpi/types.h>
#include <cedar/kernels/restrict.h>

namespace cedar { namespace cdr3 { namespace mpi {

class restrict_f90 : public kernels::restriction<stypes>
{
	void run(const restrict_op & R,
	         const grid_func & x,
	         grid_func & y) override;
};

}}}

#endif
