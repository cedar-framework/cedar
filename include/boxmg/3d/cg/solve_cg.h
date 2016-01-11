#ifndef BOXMG_3D_SOLVE_CG_H
#define BOXMG_3D_SOLVE_CG_H

#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/mpi/grid_func.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg3::mpi;
	void fortran_solve_cg(grid_func & x,
	                      const grid_func & b,
	                      const grid_func & ABD,
	                      real_t *bbd);

	void mpi_solve_cg_lu(mpi::grid_func &x,
	                     const mpi::grid_func &b,
	                     const mpi::grid_func & ABD,
	                     real_t *bbd);
}

}}}

#endif
