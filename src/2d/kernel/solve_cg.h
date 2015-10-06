#ifndef BOXMG_2D_KERNEL_SOLVE_CG_H
#define BOXMG_2D_KERNEL_SOLVE_CG_H

#include "core/grid_func.h"
#include "solver/boxmg.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void fortran_solve_cg(core::GridFunc & x,
	                      const core::GridFunc & b,
	                      const core::GridFunc & ABD,
	                      real_t *bbd);


	void mpi_solve_cg_lu(core::GridFunc &x,
	                     const core::GridFunc &b,
	                     const core::GridFunc & ABD,
	                     real_t *bbd);


	void solve_cg_boxmg(const solver::BoxMG & cg_solver,
	                    core::GridFunc &x,
	                    const core::GridFunc &b);

}

}}}

#endif
