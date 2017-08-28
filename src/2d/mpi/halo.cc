#include <cstdlib>
#include <iostream>

#include <cedar/types.h>
#include <cedar/util/timer.h>
#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"
#include "cedar/2d/ftn/mpi/BMG_parameters_c.h"
#include "cedar/2d/mpi/halo.h"

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_SETUP_MSG(int *pMSG, int *pMSGSO, len_t *iMSG_Geom,
	                           len_t NMSGi, int *pSI_MSG, int IBC, len_t *IGRD,
	                           int nog, int nogm, int nproc, int myproc,
	                           len_t* dimx, len_t *dimy, len_t *dimxfine, len_t *dimyfine,
	                           int *procgrid, int nproci, int nprocj, int mpicomm);
	void BMG2_SymStd_SETUP_LS(len_t *iWorkMSG, len_t NMSGi, int *pMSG, int *pLS, int *pSI_MSG,
	                          int *procgrid, int nproci, int nprocj, int nog);
	void BMG2_SymStd_SETUP_fine_stencil(int kf, real_t *SO, len_t IIF, len_t JJF, int NStncl,
	                                    len_t *iWork, len_t NMSGi, int *pMSGSO,
	                                    real_t *msg_buffer, len_t NMSGr, int mpicomm);
	void BMG2_SymStd_UTILS_update_ghosts(int K, real_t *x, len_t Nx, len_t Ny, len_t *iWork,
	                                     len_t NMSGi, int *pMSG,
	                                     real_t *buffer, len_t NMSGr, int nog, int mpicomm);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr2 { namespace kernel {
namespace impls
{


	std::unique_ptr<halo_exchanger> setup_msg(const kernel_params & params, grid_topo & topo)
	{
		auto halof = std::make_unique<msg_exchanger>(topo);
		auto & ctx = halof->context();
		int rank;
		int ibc;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		MPI_Comm_rank(topo.comm, &rank);
		rank++; // Fortran likes to be difficult...

		BMG_get_bc(params.per_mask(), &ibc);

		BMG2_SymStd_SETUP_MSG(ctx.pMSG.data(), ctx.pMSGSO.data(),
		                      ctx.msg_geom.data(), ctx.msg_geom.size(),
		                      &ctx.pSI_MSG, ibc, topo.IGRD(), topo.nlevel(),
		                      topo.nlevel(), topo.nproc(), rank,
		                      ctx.dimx.data(), ctx.dimy.data(),
		                      ctx.dimxfine.data(), ctx.dimyfine.data(),
		                      ctx.proc_grid.data(), topo.nproc(0), topo.nproc(1),
		                      fcomm);

		BMG2_SymStd_SETUP_LS(ctx.msg_geom.data(), ctx.msg_geom.size(),
		                     ctx.pMSG.data(), ctx.pLS.data(), &ctx.pSI_MSG,
		                     ctx.proc_grid.data(), topo.nproc(0), topo.nproc(1),
		                     topo.nlevel());

		halof->init();

		return halof;
	}


	template<>
	void msg_stencil_exchange(const kernel_params & params, halo_exchanger *halof, mpi::stencil_op<five_pt> & sop)
	{
		MsgCtx *ctx = (MsgCtx*) halof->context_ptr();
		grid_topo &topo = sop.grid();
		int nstencil;

		nstencil = 3;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		timer_begin("halo-stencil");
		BMG2_SymStd_SETUP_fine_stencil(topo.level()+1, sop.data(),
		                               sop.len(0), sop.len(1), nstencil,
		                               ctx->msg_geom.data(), ctx->msg_geom.size(),
		                               ctx->pMSGSO.data(), ctx->msg_buffer.data(),
		                               ctx->msg_buffer.size(), fcomm);
		timer_end("halo-stencil");
	}


	template<>
	void msg_stencil_exchange(const kernel_params & params, halo_exchanger *halof, mpi::stencil_op<nine_pt> & sop)
	{
		MsgCtx *ctx = (MsgCtx*) halof->context_ptr();
		grid_topo &topo = sop.grid();
		int nstencil;

		nstencil = 5;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		timer_begin("halo-stencil");
		BMG2_SymStd_SETUP_fine_stencil(topo.level()+1, sop.data(),
		                               sop.len(0), sop.len(1), nstencil,
		                               ctx->msg_geom.data(), ctx->msg_geom.size(),
		                               ctx->pMSGSO.data(), ctx->msg_buffer.data(),
		                               ctx->msg_buffer.size(), fcomm);
		timer_end("halo-stencil");
	}


	void msg_exchange(const kernel_params & params, halo_exchanger *halof, mpi::grid_func & f)
	{
		MsgCtx *ctx = (MsgCtx*) halof->context_ptr();
		grid_topo &topo = f.grid();

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		timer_begin("halo");
		BMG2_SymStd_UTILS_update_ghosts(topo.level()+1, f.data(), f.len(0), f.len(1), ctx->msg_geom.data(),
		                                ctx->msg_geom.size(), ctx->pMSG.data(), ctx->msg_buffer.data(),
		                                ctx->msg_buffer.size(), topo.nlevel(), fcomm);
		timer_end("halo");
	}
}
}}}
