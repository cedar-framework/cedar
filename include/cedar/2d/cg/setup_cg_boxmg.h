#ifndef CEDAR_2D_KERNEL_SETUP_CG_CEDAR_H
#define CEDAR_2D_KERNEL_SETUP_CG_CEDAR_H

#include <memory>

#include <cedar/kernel_params.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/solver.h>


namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr2::mpi;
	template <class sten>
	void setup_cg_boxmg(const kernel_params & params,
	                    const mpi::stencil_op<sten> & so,
	                    std::shared_ptr<config::reader> conf,
	                    std::shared_ptr<solver<sten>> *slv);
}

}}}

#endif
