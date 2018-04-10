#ifndef CEDAR_3D_INTERFACE_TYPES_H
#define CEDAR_3D_INTERFACE_TYPES_H

#include <cedar/3d/mpi/kernel_manager.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/util/topo.h>

namespace cedar { namespace cdr3 {

		struct op_container
		{
			op_container(std::shared_ptr<cedar::grid_topo> topo,
			             config::reader & conf) : op(topo), xgf(topo), bgf(topo)
			{
				kman = mpi::build_kernel_manager(conf);
			}
			mpi::stencil_op<xxvii_pt> op;
			mpi::grid_func xgf;
			mpi::grid_func bgf;
			mpi::kman_ptr kman;
		};
}}

#endif
