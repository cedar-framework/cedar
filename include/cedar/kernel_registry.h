#ifndef CEDAR_KERNEL_REGISTRY_H
#define CEDAR_KERNEL_REGISTRY_H

#include <type_traits>

#include <cedar/types.h>
#include <cedar/type_list.h>
#include <cedar/kernels/point_relax.h>
#include <cedar/kernels/line_relax.h>
#include <cedar/kernels/coarsen_op.h>
#include <cedar/kernels/interp_add.h>
#include <cedar/kernels/matvec.h>
#include <cedar/kernels/residual.h>
#include <cedar/kernels/restrict.h>
#include <cedar/kernels/setup_interp.h>
#include <cedar/kernels/solve_cg.h>

#include <cedar/kernels/halo_exchange.h>
#include <cedar/kernels/setup_nog.h>

namespace cedar
{
	template<class solver_types>
	using ser_kernels = type_list<kernels::point_relax<solver_types>,
	                              kernels::line_relax<solver_types, relax_dir::x>,
	                              kernels::line_relax<solver_types, relax_dir::y>,
	                              kernels::coarsen_op<solver_types>,
	                              kernels::restriction<solver_types>,
	                              kernels::interp_add<solver_types>,
	                              kernels::residual<solver_types>,
	                              kernels::setup_interp<solver_types>,
	                              kernels::solve_cg<solver_types>>;

	template<class solver_types>
	using mpi_kernels = type_list<kernels::halo_exchange<solver_types>,
	                              kernels::matvec<solver_types>,
	                              kernels::setup_nog<solver_types>>;

	template<class solver_types, exec_mode mode>
	using klist = typename std::conditional<mode==exec_mode::serial,
	                                        ser_kernels<solver_types>,
	                                        typename type_cat<ser_kernels<solver_types>,
	                                                          mpi_kernels<solver_types>>::type
	                                        >::type;
}

#endif
