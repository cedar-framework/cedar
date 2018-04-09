#ifndef CEDAR_3D_MPI_TYPES_H
#define CEDAR_3D_MPI_TYPES_H

#include <cedar/solver_types.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/grid_func.h>
#include <cedar/3d/mpi/restrict_op.h>
#include <cedar/3d/mpi/prolong_op.h>
#include <cedar/3d/relax_stencil.h>

namespace cedar
{
	namespace cdr3
	{
		namespace mpi
		{
			using stypes = solver_types<
				stencil_op,
				seven_pt,
				xxvii_pt,
				grid_func,
				prolong_op,
				restrict_op,
				relax_stencil>;
		}
	}
}

#endif
