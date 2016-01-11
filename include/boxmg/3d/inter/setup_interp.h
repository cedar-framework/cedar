#ifndef BOXMG_3D_INTER_SETUP_INTERP_H
#define BOXMG_3D_INTER_SETUP_INTERP_H

#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/inter/prolong_op.h>
#include <boxmg/3d/mpi/stencil_op.h>
#include <boxmg/3d/inter/mpi/prolong_op.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg3::mpi;
	void setup_interp(int kf, int kc, int nog,
	                  const stencil_op & fop,
	                  const stencil_op & cop,
	                  inter::prolong_op & P);

	void mpi_setup_interp(int kf, int kc, int nog,
	                      const mpi::stencil_op & fop,
	                      const mpi::stencil_op & cop,
	                      inter::mpi::prolong_op & P);
}

}}}

#endif
