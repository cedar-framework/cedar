#include "cedar/2d/matvec.h"
#include "cedar/2d/mpi/stencil_op.h"

#include "cedar/2d/mpi/halo.h"

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_UTILS_matvec(int k, real_t *SO, real_t *QF,
	                              real_t *Q, len_t II, len_t JJ,
	                              int kf, int ifd, int nstencil,
	                              len_t *iWork, int *pMSG,
	                              real_t *msg_buffer, int MPICOMM);
}


namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr2::mpi;
	void matvec(const mpi::stencil_op & so,
	            const mpi::grid_func & x,
	            mpi::grid_func & b)
	{
		using namespace cedar::cdr2;
		int k, kf, ifd;
		int nstencil;

		mpi::stencil_op & sod = const_cast<mpi::stencil_op&>(so);
		mpi::grid_func & xd = const_cast<mpi::grid_func&>(x);
		grid_topo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) sod.halo_ctx;
		grid_stencil & sten = sod.stencil();

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

