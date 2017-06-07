#ifndef CEDAR_3D_SETUP_CG_LU_MPI_H
#define CEDAR_3D_SETUP_CG_LU_MPI_H

#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"
#include "cedar/2d/ftn/BMG_parameters_c.h"
#include <cedar/kernel_params.h>
#include <cedar/3d/mpi/stencil_op.h>


extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_SETUP_cg_LU(real_t *so, len_t ii, len_t jj, len_t kk,
	                                 int nstencil, real_t * abd, len_t nabd1, len_t nabd2,
	                                 real_t *ws, len_t nmsgr,
	                                 int nproci, int nprocj, int nprock, int nproc, int myproc,
	                                 int *proc_grid, int *proc_coord, len_t *loc_arr_size, int mpicomm);
}

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;

	template<class sten>
	void mpi_setup_cg_lu(const kernel_params & params,
	                     const mpi::stencil_op<sten> & so,
	                     grid_func & ABD)
	{
		auto & copd = const_cast<mpi::stencil_op<sten>&>(so);

		grid_topo & topo = copd.grid();
		MsgCtx *ctx = (MsgCtx*) copd.halo_ctx;
		int nstencil = stencil_ndirs<sten>::value;

		int rank;
		MPI_Comm_rank(topo.comm, &rank);
		rank++; // 1 based indexing

		len_t local_arr_ptr = ctx->pMSG(ipL_MSG_LocalArraySize,0) - 1;  // 1 vs 0 based indexing

		if (params.per_mask()) {
			log::error << "MPI LU cg solver does not support periodic BCs" << std::endl;
		}

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_BMG3_SymStd_SETUP_cg_LU(copd.data(), copd.len(0), copd.len(1), copd.len(2),
		                            nstencil, ABD.data(), ABD.len(0), ABD.len(1),
		                            ctx->msg_buffer.data(), ctx->msg_buffer.size(),
		                            topo.nproc(0), topo.nproc(1), topo.nproc(2), topo.nproc(),
		                            rank,
		                            ctx->proc_grid.data(), ctx->proc_coord.data(),
		                            &ctx->msg_geom.data()[local_arr_ptr],
		                            fcomm);
	}
}

}}}

#endif
