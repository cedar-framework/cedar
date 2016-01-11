#ifndef BOXMG_2D_KERNEL_SETUP_RELAX_H
#define BOXMG_2D_KERNEL_SETUP_RELAX_H

#include "boxmg/2d/stencil_op.h"
#include "boxmg/2d/relax_stencil.h"
#include "boxmg/2d/mpi/stencil_op.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg2d::mpi;
	void setup_rbgs_point(const stencil_op & so,
	                      relax_stencil & sor);
	void setup_rbgs_x(const stencil_op & so,
	                  relax_stencil & sor);
	void setup_rbgs_y(const stencil_op & so,
	                  relax_stencil & sor);
	void mpi_setup_rbgs_point(const mpi::stencil_op & so,
	                          relax_stencil & sor);
	void mpi_setup_rbgs_x(const mpi::stencil_op & so,
	                      relax_stencil & sor);
	void mpi_setup_rbgs_y(const mpi::stencil_op & so,
	                      relax_stencil & sor);
}

}}}

#endif
