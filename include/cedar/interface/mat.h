#ifndef CEDAR_INTERFACE_MAT_H
#define CEDAR_INTERFACE_MAT_H

#include <memory>

#include <cedar/capi.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/kernel_manager.h>
#ifdef ENABLE_3D
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/kernel_manager.h>
#endif

struct cedar_mat_cont
{
	unsigned short nd;
	bool compressed;  /** Compressed stencil in 2D means five_pt */
	cedar_topo topo;
	std::unique_ptr<cedar::cdr2::mpi::stencil_op<cedar::cdr2::five_pt>> op2comp;
	std::unique_ptr<cedar::cdr2::mpi::stencil_op<cedar::cdr2::nine_pt>> op2full;
	cedar::cdr2::mpi::kman_ptr kman2;
	#ifdef ENABLE_3D
	std::unique_ptr<cedar::cdr3::mpi::stencil_op<cedar::cdr3::seven_pt>> op3comp;
	std::unique_ptr<cedar::cdr3::mpi::stencil_op<cedar::cdr3::xxvii_pt>> op3full;
	cedar::cdr3::mpi::kman_ptr kman3;
	#endif
};

cedar_mat_cont * cedar_mat_getobj(cedar_mat mat);

#endif
