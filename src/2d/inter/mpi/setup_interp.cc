#include "boxmg/2d/mpi/halo.h"
#include "boxmg/2d/inter/setup_interp.h"

extern "C" {
	using namespace boxmg;
	void MPI_BMG2_SymStd_SETUP_interp_OI(int kf, int kc, real_t *so, real_t *ci,
	                                     len_t iif, len_t jjf, len_t iic, len_t jjc,
	                                     int nog, int ifd, int nstncl, int nogm,
	                                     len_t *igrd, len_t *iWork, len_t NMSGi, int *pMSG,
	                                     real_t *msg_buffer, len_t nmsgr, int mpicomm);
}

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	using namespace boxmg::bmg2d;

	void mpi_setup_interp(int kf, int kc, int nog,
	                      const mpi::stencil_op & fop,
	                      const mpi::stencil_op & cop,
	                      inter::mpi::prolong_op & P)
	{
		int ifd, nstencil;

		mpi::stencil_op & fopd = const_cast<mpi::stencil_op&>(fop);
		mpi::stencil_op & copd = const_cast<mpi::stencil_op&>(cop);
		grid_stencil & fsten = fopd.stencil();
		grid_stencil & csten = copd.stencil();
		grid_topo & topo = fopd.grid();
		MsgCtx *ctx = (MsgCtx*) fopd.halo_ctx;

		if (fsten.five_pt()) {
			ifd = 1;
			nstencil = 3;
		} else {
			ifd = 0;
			nstencil = 5;
		}
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		MPI_BMG2_SymStd_SETUP_interp_OI(kf, kc, fopd.data(), P.data(),
		                                fsten.len(0), fsten.len(1), csten.len(0), csten.len(1),
		                                nog, ifd, nstencil, nog, topo.IGRD(),
		                                ctx->msg_geom.data(), ctx->msg_geom.size(),
		                                ctx->pMSG.data(), ctx->msg_buffer.data(),
		                                ctx->msg_buffer.size(), fcomm);
	}
}

}}}
