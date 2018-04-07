#ifndef CEDAR_2D_KERNEL_MPI_SETUP_CG_LU_H
#define CEDAR_2D_KERNEL_MPI_SETUP_CG_LU_H

#include <cedar/kernel_params.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/mpi/msg_exchanger.h>
#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_cg_LU(real_t* SO, len_t II, len_t JJ, int nstncl,
	                                 len_t iGs, len_t jGs, len_t nGx, len_t nGy,
	                                 real_t *ABD, len_t NABD1, len_t NABD2,
	                                 real_t *ws, len_t NMSGr, int nproci, int nprocj, int nproc,
	                                 int myproc, int *proc_grid, int *proc_coord, len_t *locarrsize,
	                                 int MPICOMM);
}

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr2::mpi;
	template <class sten>
		void mpi_setup_cg_lu(const kernel_params & params, mpi::msg_exchanger * halof,
		                     const mpi::stencil_op<sten> & so,
		                     grid_func & ABD)
	{
		auto & copd = const_cast<mpi::stencil_op<sten>&>(so);

		grid_topo & topo = copd.grid();
		MsgCtx *ctx = (MsgCtx*) halof->context_ptr();
		int nstencil;

		nstencil = stencil_ndirs<sten>::value;

		int rank;
		MPI_Comm_rank(topo.comm, &rank);
		rank++; // 1 based indexing

		len_t local_arr_ptr = ctx->pMSG(ipL_MSG_LocalArraySize,0) - 1;  // 1 vs 0 based indexing

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		if (params.per_mask()) {
			log::error << "MPI LU cg solver does not support periodic BCs" << std::endl;
		}

		MPI_BMG2_SymStd_SETUP_cg_LU(copd.data(), so.len(0), so.len(1),
		                            nstencil, topo.is(0), topo.is(1),
		                            topo.nglobal(0), topo.nglobal(1),
		                            ABD.data(), ABD.len(0), ABD.len(1),
		                            ctx->msg_buffer.data(),
		                            ctx->msg_buffer.size(),
		                            topo.nproc(0), topo.nproc(1), topo.nproc(),
		                            rank,
		                            ctx->proc_grid.data(),
		                            ctx->proc_coord.data(),
		                            &ctx->msg_geom.data()[local_arr_ptr],
		                            fcomm);

	}
}

}}}

#endif
