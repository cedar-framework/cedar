#include "halo.h"
#include "setup_interp.h"

extern "C" {
	void bmg2_symstd_setup_interp_oi(int*, int*, double*, double*, double*,
	                                  int*, int*, int*, int*, int*, int*, int*, int*);

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

	void setup_interp(int kf, int kc, int nog,
	                  const StencilOp &fop,
	                  const StencilOp &cop,
	                  inter::prolong_op &P)
	{
		using namespace boxmg::bmg2d;
		int ifd;
		int iif, jjf, iic, jjc;
		int nstencil;
		int irelax = 0;

		const GridStencil &fsten = fop.stencil();
		const GridStencil &csten = cop.stencil();
		StencilOp &fopd = const_cast<StencilOp&>(fop);
		StencilOp &copd = const_cast<StencilOp&>(cop);

		P.fine_op = &fopd;

		iif = fsten.len(0);
		jjf = fsten.len(1);
		iic = csten.len(0);
		jjc = csten.len(1);

		if (fsten.five_pt()) {
			ifd = 0;
			nstencil = 3;
		} else {
			ifd = 1;
			nstencil = 5;
		}

		// std::cout << "fine stencil: " << fopd.data() << std::endl;
		// std::cout << "coarse stencil: " << copd.data() << std::endl;
		// std::cout << "interpolation op: " << P.data() << std::endl;
		bmg2_symstd_setup_interp_oi(&kf, &kc, fopd.data(), copd.data(), P.data(), &iif, &jjf, &iic, &jjc,
		                            &nog, &ifd, &nstencil, &irelax);
	}

	void mpi_setup_interp(int kf, int kc, int nog,
	                      const StencilOp & fop,
	                      const StencilOp & cop,
	                      inter::prolong_op & P)
	{
		int ifd, nstencil;

		mpi::StencilOp & fopd = const_cast<mpi::StencilOp&>(dynamic_cast<const mpi::StencilOp&>(fop));
		mpi::StencilOp & copd = const_cast<mpi::StencilOp&>(dynamic_cast<const mpi::StencilOp&>(cop));
		GridStencil & fsten = fopd.stencil();
		GridStencil & csten = copd.stencil();
		mpi::GridTopo & topo = fopd.grid();
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
