#ifndef BOXMG_2D_KERNEL_RELAX_H
#define BOXMG_2D_KERNEL_RELAX_H

#include "core/stencil_op.h"
#include "core/relax_stencil.h"
#include "core/grid_func.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void relax_rbgs_point(const StencilOp & so,
	                      grid_func & x,
	                      const grid_func & b,
	                      const RelaxStencil & sor,
	                      cycle::Dir cycle_dir);

	void relax_lines_x(const StencilOp & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const RelaxStencil & sor,
	                   grid_func & res,
	                   cycle::Dir cycle_dir);

	void relax_lines_y(const StencilOp & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const RelaxStencil & sor,
	                   grid_func & res,
	                   cycle::Dir cycle_dir);

	void mpi_relax_rbgs_point(const StencilOp & so,
	                          grid_func & x,
	                          const grid_func & b,
	                          const RelaxStencil & sor,
	                          cycle::Dir cycle_dir);

	void mpi_relax_lines_x(const StencilOp & so,
	                       grid_func & x,
	                       const grid_func & b,
	                       const RelaxStencil & sor,
	                       grid_func & res,
	                       cycle::Dir cycle_dir);

	void mpi_relax_lines_y(const StencilOp & so,
	                       grid_func & x,
	                       const grid_func & b,
	                       const RelaxStencil & sor,
	                       grid_func & res,
	                       cycle::Dir cycle_dir);
}

}}}

#endif
