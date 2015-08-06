#ifndef BOXMG_2D_KERNEL_SETUP_RELAX_H
#define BOXMG_2D_KERNEL_SETUP_RELAX_H

#include "core/stencil_op.h"
#include "core/relax_stencil.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void setup_rbgs_point(const core::StencilOp & so,
	                      core::RelaxStencil & sor);
	void mpi_setup_rbgs_point(const core::StencilOp & so,
	                          core::RelaxStencil & sor);
}

}}}

#endif
