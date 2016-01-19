#ifndef BOXMG_3D_GALERKIN_PROD_H
#define BOXMG_3D_GALERKIN_PROD_H

#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/inter/prolong_op.h>
#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/mpi/stencil_op.h>
#include <boxmg/3d/inter/mpi/prolong_op.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg3::mpi;
	void galerkin_prod(int kf, int kc, int nog,
	                   const inter::prolong_op &P,
	                   const stencil_op & fop,
	                   stencil_op & cop);

	void mpi_galerkin_prod(int kf, int kc, int nog,
	                       const inter::mpi::prolong_op & P,
	                       const mpi::stencil_op & fop,
	                       mpi::stencil_op & cop);
}

}}}

#endif