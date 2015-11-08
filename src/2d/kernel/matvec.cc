#include "matvec.h"
#include "core/mpi/stencil_op.h"

#include "halo.h"

extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_UTILS_matvec(int k, real_t *SO, real_t *QF,
	                              real_t *Q, len_t II, len_t JJ,
	                              int kf, int ifd, int nstencil,
	                              len_t *iWork, int *pMSG,
	                              real_t *msg_buffer, int MPICOMM);
}


namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void matvec(const StencilOp & so,
	            const grid_func & x,
	            grid_func & b)
	{
		using namespace boxmg::bmg2d;
		int k, kf, ifd;
		int nstencil;

		mpi::StencilOp & sod = const_cast<mpi::StencilOp&>(dynamic_cast<const mpi::StencilOp&>(so));
		grid_func & xd = const_cast<grid_func&>(x);
		mpi::GridTopo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) sod.halo_ctx;
		GridStencil & sten = sod.stencil();

		k = topo.level()+1;
		kf = topo.nlevel();
		if (sten.five_pt()) {
			ifd = 1;
			nstencil = 3;
		} else {
			ifd = 0;
			nstencil = 5;
		}

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		BMG2_SymStd_UTILS_matvec(k, sod.data(), b.data(),
		                         xd.data(), sten.len(0),
		                         sten.len(1), kf, ifd, nstencil,
		                         ctx->msg_geom.data(),
		                         ctx->pMSG.data(),
		                         ctx->msg_buffer.data(),
		                         fcomm);
	}
}

}}}

