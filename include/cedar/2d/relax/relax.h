#ifndef CEDAR_2D_KERNEL_RELAX_H
#define CEDAR_2D_KERNEL_RELAX_H

#include <cedar/kernel_params.h>
#include <cedar/halo_exchanger.h>
#include "cedar/cycle/types.h"
#include "cedar/2d/stencil_op.h"
#include "cedar/2d/relax_stencil.h"
#include "cedar/2d/grid_func.h"

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	template<class sten>
	void relax_rbgs_point(const kernel_params & params,
	                      const stencil_op<sten> & so,
	                      grid_func & x,
	                      const grid_func & b,
	                      const relax_stencil & sor,
	                      cycle::Dir cycle_dir);

	template<class sten>
	void relax_lines_x(const kernel_params & params,
	                   const stencil_op<sten> & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func & res,
	                   cycle::Dir cycle_dir);

	template<class sten>
	void relax_lines_y(const kernel_params & params,
	                   const stencil_op<sten> & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func & res,
	                   cycle::Dir cycle_dir);
}

}}}

#endif
