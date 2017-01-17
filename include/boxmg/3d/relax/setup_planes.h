#ifndef BOXMG_3D_SETUP_PLANES_H
#define BOXMG_3D_SETUP_PLANES_H

#include <boxmg/2d/solver.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/relax_stencil.h>
#include <boxmg/3d/mpi/stencil_op.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	using slv2_ptr = std::unique_ptr<::boxmg::bmg2d::solver>;
	void setup_relax_xy(const stencil_op & so,
	                    std::vector<slv2_ptr> & planes);
	void setup_relax_xz(const stencil_op & so,
	                    std::vector<slv2_ptr> & planes);
	void setup_relax_yz(const stencil_op & so,
	                    std::vector<slv2_ptr> & planes);
}

}}}


#endif