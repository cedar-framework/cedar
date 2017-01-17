#ifndef BOXMG_3D_RELAX_PLANES_H
#define BOXMG_3D_RELAX_PLANES_H

#include <boxmg/2d/solver.h>
#include <boxmg/cycle/types.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/grid_func.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	using slv2_ptr = std::unique_ptr<::boxmg::bmg2d::solver>;

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
