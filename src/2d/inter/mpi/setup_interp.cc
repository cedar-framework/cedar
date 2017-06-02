#include "cedar/2d/mpi/halo.h"
#include "cedar/2d/inter/setup_interp.h"

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_interp_OI(int kf, int kc, real_t *so, real_t *ci,
	                                     len_t iif, len_t jjf, len_t iic, len_t jjc,
	                                     int nog, int ifd, int nstncl, int nogm, int ibc,
	                                     len_t *igrd, len_t *iWork, len_t NMSGi, int *pMSG,
	                                     real_t *msg_buffer, len_t nmsgr, int mpicomm);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr2;

	template<>
	void mpi_setup_interp(const kernel_params & params,
	                      const mpi::stencil_op<five_pt> & fop,
	                      const mpi::stencil_op<nine_pt> & cop,
	                      inter::mpi::prolong_op & P)
	{
		int ifd, nstencil, jpn;
		int kf, kc, nog;

		auto & fopd = const_cast<mpi::stencil_op<five_pt>&>(fop);
		grid_topo & topo = fopd.grid();
		MsgCtx *ctx = (MsgCtx*) fopd.halo_ctx;

		ifd = 1;
		nstencil = 3;

		kc = topo.level() + 1;
		nog = topo.nlevel();
		kf = kc + 1;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		BMG_get_bc(params.per_mask(), &jpn);

		MPI_BMG2_SymStd_SETUP_interp_OI(kf, kc, fopd.data(), P.data(),
		                                fop.len(0), fop.len(1), cop.len(0), cop.len(1),
		                                nog, ifd, nstencil, nog, jpn, topo.IGRD(),
		                                ctx->msg_geom.data(), ctx->msg_geom.size(),
		                                ctx->pMSG.data(), ctx->msg_buffer.data(),
		                                ctx->msg_buffer.size(), fcomm);
	}


	template<>
	void mpi_setup_interp(const kernel_params & params,
	                      const mpi::stencil_op<nine_pt> & fop,
	                      const mpi::stencil_op<nine_pt> & cop,
	                      inter::mpi::prolong_op & P)
	{
		int ifd, nstencil, jpn;
		int kf, kc, nog;

		auto & fopd = const_cast<mpi::stencil_op<nine_pt>&>(fop);

		grid_topo & topo = fopd.grid();
		MsgCtx *ctx = (MsgCtx*) fopd.halo_ctx;

		ifd = 0;
		nstencil = 5;

		kc = topo.level() + 1;
		nog = topo.nlevel();
		kf = kc + 1;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		BMG_get_bc(params.per_mask(), &jpn);

		MPI_BMG2_SymStd_SETUP_interp_OI(kf, kc, fopd.data(), P.data(),
		                                fop.len(0), fop.len(1), cop.len(0), cop.len(1),
		                                nog, ifd, nstencil, nog, jpn, topo.IGRD(),
		                                ctx->msg_geom.data(), ctx->msg_geom.size(),
		                                ctx->pMSG.data(), ctx->msg_buffer.data(),
		                                ctx->msg_buffer.size(), fcomm);
	}
}

}}}
