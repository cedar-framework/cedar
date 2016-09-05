#include <cstdlib>
#include <iostream>

#include <boxmg/types.h>
#include "boxmg/2d/ftn/mpi/BMG_workspace_c.h"
#include "boxmg/2d/ftn/mpi/BMG_parameters_c.h"
#include "boxmg/2d/mpi/halo.h"

extern "C" {
	using namespace boxmg;
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
}

namespace boxmg { namespace bmg2d { namespace kernel {
namespace impls
{
	MsgCtx::MsgCtx(grid_topo & topo) :
		pMSG(NBMG_pMSG, topo.nlevel()),
		pLS(NBMG_pLS, topo.nlevel()),
		pMSGSO(NBMG_pMSG, topo.nlevel()),
		proc_grid(topo.nproc(0), topo.nproc(1)),
		proc_coord(topo.nproc(0)*topo.nproc(1)*2),
		dimxfine(topo.nproc(0)),
		dimyfine(topo.nproc(1)),
		dimx(topo.nproc(0), topo.nlevel()),
		dimy(topo.nproc(1), topo.nlevel())
	{
		MPI_Comm_split(topo.comm, topo.coord(1), topo.coord(0),
		               &xlinecomm);
		MPI_Comm_split(topo.comm, topo.coord(0), topo.coord(1),
		               &ylinecomm);
		len_t NLx = topo.nlocal(0) - 2;
		len_t NLy = topo.nlocal(1) - 2;
		// !
		// !  Storage for NLx, NLy, NLz arrays for all levels
		// !
		len_t NMSGi = (topo.nproc(0)+topo.nproc(1))*topo.nlevel();
		// !
		// !  Add storage for MSG workspace (Not a sharp bound)
		// !
		NMSGi = NMSGi + topo.nproc()+2 + topo.nlevel()*(16*(NLx+NLy+6)+18*topo.nproc()+26);
		// !
		// ! Add storage for MSGSO workspace
		// !
		NMSGi = NMSGi + topo.nproc()+2 + topo.nlevel()*(24*(NLx+NLy+4)+18*topo.nproc()+26);
		// !
		// ! Add storage for Line Solves
		// !
		NMSGi = NMSGi + 4*topo.nlevel()*topo.nproc();
		msg_geom.resize(NMSGi);

		// !
		// !  Need to fix the this bound!!
		// !
		// !  - MSG real buffer space 
		// !  - Workspace for coarse-grid solve communication.
		// !
		// NMSGr = std::max(2*std::max(topo.nproc(0)+6, topo.nproc(1)+6),
		//                  10*(5*NGx_c*NGy_c+5),
		//               &              4*NLy*NProcI + 5*NLx, 4*NLx*NProcJ + 5*NLy )

		p_NLx_kg = 1;
		p_NLy_kg = p_NLx_kg + topo.nproc(0)*topo.nlevel();
		pSI_MSG = p_NLy_kg + topo.nproc(1)*topo.nlevel();
		grid_topo cg_topo(topo.igrd_ptr(), 0, topo.nlevel());
		len_t NGx_c = cg_topo.nglobal(0) - 2;
		len_t NGy_c = cg_topo.nglobal(1) - 2;
		len_t NMSGr = std::max(4*NLy*topo.nproc(0) + 5*NLx,
		                       4*NLx*topo.nproc(1) + 5*NLy);
		NMSGr = std::max(NMSGr, 10*(5*NGx_c*NGy_c+5));
		NMSGr = std::max(NMSGr, 2*std::max(NLx+6, NLy+6));
		msg_buffer.resize(NMSGr);

		auto pc = proc_coord.begin();
		for (auto j : range(topo.nproc(1))) {
			dimyfine[j] = topo.nlocal(1) - 2;
			for (auto i : range(topo.nproc(0))) {
				proc_grid(i,j) = j*topo.nproc(0) + i + 1;
				*pc = i+1;
				pc++;
				*pc = j+1;
				pc++;
			}
		}

		if (topo.dimyfine.size() > 0) {
			dimyfine = topo.dimyfine;
		} else {
			for (auto j : range(topo.nproc(1))) {
				dimyfine[j] = topo.nlocal(1) - 2;
			}
		}

		if (topo.dimxfine.size() > 0) {
			dimxfine = topo.dimxfine;
		} else {
			for (auto i : range(topo.nproc(0))) {
				dimxfine[i] = topo.nlocal(0) - 2;
			}
		}
	}


	void setup_msg(grid_topo & topo, void **msg_ctx)
	{
		MsgCtx *ctx = new MsgCtx(topo);
		int rank;
		int ibc;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		MPI_Comm_rank(topo.comm, &rank);
		rank++; // Fortran likes to be difficult...

		ibc = BMG_BCs_definite;

		BMG2_SymStd_SETUP_MSG(ctx->pMSG.data(), ctx->pMSGSO.data(),
		                      ctx->msg_geom.data(), ctx->msg_geom.size(),
		                      &ctx->pSI_MSG, ibc, topo.IGRD(), topo.nlevel(),
		                      topo.nlevel(), topo.nproc(), rank,
		                      ctx->dimx.data(), ctx->dimy.data(),
		                      ctx->dimxfine.data(), ctx->dimyfine.data(),
		                      ctx->proc_grid.data(), topo.nproc(0), topo.nproc(1),
		                      fcomm);

		BMG2_SymStd_SETUP_LS(ctx->msg_geom.data(), ctx->msg_geom.size(),
		                     ctx->pMSG.data(), ctx->pLS.data(), &ctx->pSI_MSG,
		                     ctx->proc_grid.data(), topo.nproc(0), topo.nproc(1),
		                     topo.nlevel());

		*msg_ctx = (void*) ctx;
	}


	void msg_stencil_exchange(mpi::stencil_op & sop)
	{
		MsgCtx *ctx = (MsgCtx*) sop.halo_ctx;
		grid_topo &topo = sop.grid();
		grid_stencil & sten = sop.stencil();
		int nstencil;

		if (sten.five_pt()) {
			nstencil = 3;
		} else {
			nstencil = 5;
		}

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		timer_begin("halo-stencil");
		BMG2_SymStd_SETUP_fine_stencil(topo.level()+1, sop.data(),
		                               sten.len(0), sten.len(1), nstencil,
		                               ctx->msg_geom.data(), ctx->msg_geom.size(),
		                               ctx->pMSGSO.data(), ctx->msg_buffer.data(),
		                               ctx->msg_buffer.size(), fcomm);
		timer_end("halo-stencil");
	}


	void msg_exchange(mpi::grid_func & f)
	{
		MsgCtx *ctx = (MsgCtx*) f.halo_ctx;
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
