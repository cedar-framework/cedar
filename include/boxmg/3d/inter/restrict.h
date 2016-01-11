#ifndef BOXMG_3D_RESTRICT_H
#define BOXMG_3D_RESTRICT_H

#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/inter/prolong_op.h>
#include <boxmg/3d/inter/restrict_op.h>
#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/inter/mpi/restrict_op.h>
#include <boxmg/3d/mpi/grid_func.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg3::mpi;

	void fortran_restrict(const inter::restrict_op & so,
	                      const grid_func & fine,
	                      grid_func & coarse);
	void mpi_fortran_restrict(const inter::mpi::restrict_op & so,
	                          const mpi::grid_func & fine,
	                          mpi::grid_func & coarse);
}

}}}

#endif
