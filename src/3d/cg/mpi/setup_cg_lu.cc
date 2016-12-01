#include "boxmg/2d/ftn/mpi/BMG_workspace_c.h"
#include "boxmg/2d/ftn/BMG_parameters_c.h"
#include <boxmg/3d/mpi/halo.h>

#include <boxmg/3d/cg/setup_cg_lu.h>

extern "C" {
	using namespace boxmg;
	void MPI_BMG3_SymStd_SETUP_cg_LU(real_t *so, len_t ii, len_t jj, len_t kk,
	                                 int nstencil, real_t * abd, len_t nabd1, len_t nabd2,
	                                 real_t *ws, len_t nmsgr,
	                                 int nproci, int nprocj, int nprock, int nproc, int myproc,
	                                 int *proc_grid, int *proc_coord, len_t *loc_arr_size, int mpicomm);
}


namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void mpi_setup_cg_lu(const mpi::stencil_op & so,
	                     grid_func & ABD)
	{
		mpi::stencil_op & copd = const_cast<mpi::stencil_op&>(so);

		grid_stencil & csten = copd.stencil();
		grid_topo & topo = copd.grid();
		MsgCtx *ctx = (MsgCtx*) copd.halo_ctx;
		int nstencil;

		if (csten.five_pt()) {
			nstencil = 4;
		} else {
			nstencil = 14;
		}

		int rank;
		MPI_Comm_rank(topo.comm, &rank);
		rank++; // 1 based indexing

		len_t local_arr_ptr = ctx->pMSG(ipL_MSG_LocalArraySize,0) - 1;  // 1 vs 0 based indexing

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_BMG3_SymStd_SETUP_cg_LU(csten.data(), csten.len(0), csten.len(1), csten.len(2),
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
