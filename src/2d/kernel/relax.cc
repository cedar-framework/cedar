#include "fortran/BMG_parameters_c.h"
#include "core/mpi/stencil_op.h"

#include "halo.h"
#include "relax.h"

extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_relax_GS(int, real_t*, real_t*, real_t*, real_t*, len_t, len_t,
	                          int, int, int, int, int, int, int);
	void MPI_BMG2_SymStd_relax_GS(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                              len_t II, len_t JJ, int kf, int ifd, int nstncl, int irelax_sym,
	                              int updown, len_t iGs, len_t jGs, len_t *iWork, len_t NMSGi,
	                              int *pMSG, real_t *msg_buffer, len_t NMSGr, int MPICOMM);
}

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void relax_rbgs_point(const core::StencilOp & so,
	                      core::GridFunc & x,
	                      const core::GridFunc & b,
	                      const core::RelaxStencil & sor,
	                      cycle::Dir cycle_dir)
	{
		using namespace boxmg::bmg2d::core;
		int k, kf, ifd;
		int updown, nsorv, ibc, nstencil;

		const GridStencil &so_sten = so.stencil();
		StencilOp & sod = const_cast<StencilOp&>(so);
		core::RelaxStencil & sord = const_cast<core::RelaxStencil&>(sor);
		core::GridFunc & bd = const_cast<core::GridFunc&>(b);

		k = kf = 1;
		if (so_sten.five_pt()) {
			ifd = 1;
			nstencil = 3;
		} else {
			ifd = 0;
			nstencil = 5;
		}

		nsorv = 2;

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		ibc = BMG_BCs_definite;

		BMG2_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                     so_sten.len(0), so_sten.len(1), kf, ifd, nstencil, nsorv,
		                     BMG_RELAX_SYM, updown, ibc);
	}


	void mpi_relax_rbgs_point(const core::StencilOp & so,
	                          core::GridFunc & x,
	                          const core::GridFunc & b,
	                          const core::RelaxStencil & sor,
	                          cycle::Dir cycle_dir)
	{
		using namespace boxmg::bmg2d::core;
		int k, kf, ifd;
		int updown, nstencil;

		core::mpi::StencilOp & sod = const_cast<core::mpi::StencilOp&>(dynamic_cast<const core::mpi::StencilOp&>(so));
		core::mpi::GridTopo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) sod.halo_ctx;
		GridStencil & sten = sod.stencil();
		core::RelaxStencil & sord = const_cast<core::RelaxStencil&>(sor);
		core::GridFunc & bd = const_cast<core::GridFunc&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();
		if (sten.five_pt()) {
			ifd = 1;
			nstencil = 3;
		} else {
			ifd = 0;
			nstencil = 5;
		}

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		// ibc = BMG_BCs_definite;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		MPI_BMG2_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                         sten.len(0), sten.len(1), kf, ifd, nstencil, BMG_RELAX_SYM,
		                         updown, topo.is(0), topo.is(1), ctx->msg_geom.data(),
		                         ctx->msg_geom.size(), ctx->pMSG.data(), ctx->msg_buffer.data(),
		                         ctx->msg_buffer.size(), fcomm);
	}
}

}}}
