#ifndef CEDAR_2D_INTERFACE_TYPES_H
#define CEDAR_2D_INTERFACE_TYPES_H


#include <cedar/config.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/kernel_manager.h>


namespace cedar { namespace cdr2 {

		struct op_container
		{
			op_container(std::shared_ptr<cedar::grid_topo> topo,
			             config & conf) : op(topo), xgf(topo), bgf(topo)
			{
				kman = mpi::build_kernel_manager(conf);
			}
			mpi::stencil_op<nine_pt> op;
			mpi::grid_func xgf;
			mpi::grid_func bgf;
			mpi::kman_ptr kman;
		};

}}

#endif
