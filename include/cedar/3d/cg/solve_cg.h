#ifndef CEDAR_3D_SOLVE_CG_H
#define CEDAR_3D_SOLVE_CG_H

#include <cedar/3d/grid_func.h>
#include <cedar/3d/mpi/grid_func.h>
#include "cedar/3d/mpi/redist_solver.h"

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;
	void fortran_solve_cg(grid_func & x,
	                      const grid_func & b,
	                      const grid_func & ABD,
	                      real_t *bbd);

	void mpi_solve_cg_lu(mpi::grid_func &x,
	                     const mpi::grid_func &b,
	                     const mpi::grid_func & ABD,
	                     real_t *bbd);

	void solve_cg_redist(const mpi::redist_solver & cg_solver,
	                     mpi::grid_func &x,
	                     const mpi::grid_func &b);
}

}}}

#endif
