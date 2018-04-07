#ifndef CEDAR_3D_MPI_INTERP_H
#define CEDAR_3D_MPI_INTERP_H

#include <cedar/kernel_params.h>
#include <cedar/halo_exchanger.h>
#include <cedar/3d/mpi/grid_func.h>
#include <cedar/3d/inter/mpi/prolong_op.h>


namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;

	void mpi_fortran_interp(const kernel_params & params,
	                        halo_exchanger_base *halof,
	                        const inter::mpi::prolong_op & P,
	                        const mpi::grid_func & coarse,
	                        const mpi::grid_func & residual,
	                        mpi::grid_func & fine);
}

}}}


#endif
