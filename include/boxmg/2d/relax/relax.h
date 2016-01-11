#ifndef BOXMG_2D_KERNEL_RELAX_H
#define BOXMG_2D_KERNEL_RELAX_H

#include "boxmg/cycle/types.h"
#include "boxmg/2d/stencil_op.h"
#include "boxmg/2d/relax_stencil.h"
#include "boxmg/2d/grid_func.h"
#include "boxmg/2d/mpi/grid_func.h"
#include "boxmg/2d/mpi/stencil_op.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg2d::mpi;
	void relax_rbgs_point(const stencil_op & so,
	                      grid_func & x,
	                      const grid_func & b,
	                      const relax_stencil & sor,
	                      cycle::Dir cycle_dir);

	void relax_lines_x(const stencil_op & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func & res,
	                   cycle::Dir cycle_dir);

	void relax_lines_y(const stencil_op & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func & res,
	                   cycle::Dir cycle_dir);

	void mpi_relax_rbgs_point(const mpi::stencil_op & so,
	                          mpi::grid_func & x,
	                          const mpi::grid_func & b,
	                          const relax_stencil & sor,
	                          cycle::Dir cycle_dir);

	void mpi_relax_lines_x(const mpi::stencil_op & so,
	                       mpi::grid_func & x,
	                       const mpi::grid_func & b,
	                       const relax_stencil & sor,
	                       mpi::grid_func & res,
	                       cycle::Dir cycle_dir);

	void mpi_relax_lines_y(const mpi::stencil_op & so,
	                       mpi::grid_func & x,
	                       const mpi::grid_func & b,
	                       const relax_stencil & sor,
	                       mpi::grid_func & res,
	                       cycle::Dir cycle_dir);
}

}}}

#endif
