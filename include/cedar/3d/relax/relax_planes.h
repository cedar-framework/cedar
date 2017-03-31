#ifndef CEDAR_3D_RELAX_PLANES_H
#define CEDAR_3D_RELAX_PLANES_H

#include <cedar/2d/solver.h>
#include <cedar/cycle/types.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/grid_func.h>

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	using slv2_ptr = std::unique_ptr<::cedar::cdr2::solver>;

	void relax_xy(const stencil_op & so,
	              grid_func & x,
	              const grid_func & b,
	              cycle::Dir cycle_dir,
	              std::vector<slv2_ptr> & planes);

	void relax_xz(const stencil_op & so,
	              grid_func & x,
	              const grid_func & b,
	              cycle::Dir cycle_dir,
	              std::vector<slv2_ptr> & planes);

	void relax_yz(const stencil_op & so,
	              grid_func & x,
	              const grid_func & b,
	              cycle::Dir cycle_dir,
	              std::vector<slv2_ptr> & planes);

}

}}}


#endif
