#ifndef BOXMG_2D_KERNEL_SOLVE_CG_H
#define BOXMG_2D_KERNEL_SOLVE_CG_H

#include "core/grid_func.h"
#include "core/solver.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void fortran_solve_cg(GridFunc & x,
	                      const GridFunc & b,
	                      const GridFunc & ABD,
	                      real_t *bbd);


	void mpi_solve_cg_lu(GridFunc &x,
	                     const GridFunc &b,
	                     const GridFunc & ABD,
	                     real_t *bbd);


	void solve_cg_boxmg(const solver & cg_solver,
	                    GridFunc &x,
	                    const GridFunc &b);

}

}}}

#endif
