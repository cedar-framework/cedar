#ifndef CEDAR_2D_MPI_TYPES_H
#define CEDAR_2D_MPI_TYPES_H

#include <cedar/solver_types.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/inter/mpi/restrict_op.h>
#include <cedar/2d/inter/mpi/prolong_op.h>
#include <cedar/2d/relax_stencil.h>

namespace cedar
{
	namespace cdr2
	{
		namespace mpi
		{
			using stypes = solver_types<
				stencil_op,
				grid_func,
				inter::mpi::prolong_op,
				inter::mpi::restrict_op,
				relax_stencil>;
		}
	}
}

#endif
