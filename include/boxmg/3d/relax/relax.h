#ifndef BOXMG_3D_RELAX_H
#define BOXMG_3D_RELAX_H

#include <boxmg/2d/solver.h>
#include <boxmg/cycle/types.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/relax_stencil.h>
#include <boxmg/3d/mpi/grid_func.h>
#include <boxmg/3d/mpi/stencil_op.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	using slv2_ptr = std::unique_ptr<::boxmg::bmg2d::solver>;
	namespace mpi = boxmg::bmg3::mpi;
	void relax_rbgs_point(const stencil_op & so,
	                      grid_func &x,
	                      const grid_func &b,
	                      const relax_stencil & sor,
	                      cycle::Dir cycle_dir);

	void relax_xy(const stencil_op & so,
	              grid_func & x,
	              const grid_func & b,
	              cycle::Dir cycle_dir,
	              std::vector<slv2_ptr> & planes);

	void mpi_relax_rbgs_point(const mpi::stencil_op & so,
	                          mpi::grid_func & x,
	                          const mpi::grid_func & b,
	                          const relax_stencil & sor,
	                          cycle::Dir cycle_dir);
}

}}}

#endif
