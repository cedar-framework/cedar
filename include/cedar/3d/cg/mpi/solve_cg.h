#ifndef CEDAR_3D_SOLVE_CG_MPI_H
#define CEDAR_3D_SOLVE_CG_MPI_H

#include <cedar/3d/mpi/grid_func.h>

namespace cedar { namespace cdr3 { namespace mpi {
			class redist_solver;
}}}


namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;

	void mpi_solve_cg_lu(const kernel_params & params,
	                     mpi::grid_func &x,
	                     const mpi::grid_func &b,
	                     const mpi::grid_func & ABD,
	                     real_t *bbd);

	void solve_cg_redist(const kernel_params & params,
	                     const mpi::redist_solver & cg_solver,
	                     mpi::grid_func &x,
	                     const mpi::grid_func &b);
}

}}}

#endif

