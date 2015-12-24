#ifndef BOXMG_2D_KERNEL_SOLVE_CG_H
#define BOXMG_2D_KERNEL_SOLVE_CG_H

#include "boxmg/2d/grid_func.h"
#include "boxmg/2d/mpi/grid_func.h"
#include "boxmg/2d/solver.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void fortran_solve_cg(grid_func & x,
	                      const grid_func & b,
	                      const grid_func & ABD,
	                      real_t *bbd);


	void mpi_solve_cg_lu(mpi::grid_func &x,
	                     const mpi::grid_func &b,
	                     const mpi::grid_func & ABD,
	                     real_t *bbd);


	void solve_cg_boxmg(const solver & cg_solver,
	                    mpi::grid_func &x,
	                    const mpi::grid_func &b);

}

}}}

#endif
