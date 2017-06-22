#ifndef CEDAR_2D_KERNEL_HALO_H
#define CEDAR_2D_KERNEL_HALO_H

#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"

#include <cedar/kernel_params.h>
#include <cedar/array.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/grid_func.h>

namespace cedar { namespace cdr2 { namespace kernel {
namespace impls
{
	std::unique_ptr<halo_exchanger>
		setup_msg(const kernel_params & params, grid_topo &topo, void **msg_ctx);
	void msg_exchange(const kernel_params & params, mpi::grid_func & f);
	void msg_exchange(const kernel_params & params, int k, int nog, real_t *so_data, len_t nx, len_t ny, void *msg_ctx);

	void msg_stencil_exchange(const kernel_params & params, int k, int nog, real_t *so_data, len_t nx, len_t ny, void *msg_ctx);
	template<class sten>
	void msg_stencil_exchange(const kernel_params & params, mpi::stencil_op<sten> & so);
}
}}}

#endif
