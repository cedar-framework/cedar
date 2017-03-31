#ifndef CEDAR_3D_SETUP_PLANES_H
#define CEDAR_3D_SETUP_PLANES_H

#include <cedar/2d/solver.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/relax_stencil.h>
#include <cedar/3d/mpi/stencil_op.h>

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	using slv2_ptr = std::unique_ptr<::cedar::cdr2::solver>;
	void setup_relax_xy(const stencil_op & so,
	                    std::vector<slv2_ptr> & planes);
	void setup_relax_xz(const stencil_op & so,
	                    std::vector<slv2_ptr> & planes);
	void setup_relax_yz(const stencil_op & so,
	                    std::vector<slv2_ptr> & planes);
}

}}}


#endif
