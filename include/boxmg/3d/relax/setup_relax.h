#ifndef BOXMG_3D_SETUP_RELAX_H
#define BOXMG_3D_SETUP_RELAX_H

#include <boxmg/2d/solver.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/relax_stencil.h>
#include <boxmg/3d/mpi/stencil_op.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg3::mpi;
	void setup_rbgs_point(const stencil_op & so,
	                      relax_stencil & sor);
	void mpi_setup_rbgs_point(const mpi::stencil_op & so,
	                          relax_stencil & sor);
	void setup_relax_xy(const stencil_op & so,
	                    std::vector<::boxmg::bmg2d::solver> & planes);
}

}}}

#endif
