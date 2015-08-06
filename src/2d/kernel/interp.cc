#include "fortran/BMG_parameters_c.h"
#include "inter/mpi/prolong_op.h"
#include "halo.h"

#include "interp.h"

extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_interp_add(real_t*, real_t*, real_t*,real_t*,real_t*,
	                            len_t, len_t, len_t, len_t, int, int);
	void MPI_BMG2_SymStd_interp_add(len_t kc, len_t kf, len_t nog,
	                                real_t *Q, real_t *QC, real_t *RES, real_t *SO,
	                                int nstncl, real_t *CI, len_t iic, len_t jjc,
	                                len_t iif, len_t jjf, len_t iGs, len_t jGs,
	                                len_t *iWorkMSG, len_t NMSGi, int *pMSG,
	                                real_t *msg_buffer, len_t nmsgr, int mpicomm);
}

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void fortran_interp(const inter::ProlongOp & P,
	                    const core::GridFunc & coarse,
	                    const core::GridFunc & residual,
	                    core::GridFunc & fine)
	{
		using namespace boxmg::bmg2d;

		int nstencil, ibc;

		inter::ProlongOp & Pd = const_cast<inter::ProlongOp&>(P);
		core::GridFunc & coarsed = const_cast<core::GridFunc&>(coarse);
		core::GridFunc & res = const_cast<core::GridFunc&>(residual);

		if (Pd.stencil().five_pt()) {
			nstencil = 3;
		} else {
			nstencil = 5;
		}

		ibc = BMG_BCs_definite;

		BMG2_SymStd_interp_add(fine.data(), coarsed.data(), res.data(), Pd.fine_op->data(), Pd.data(),
		                       coarsed.len(0), coarsed.len(1), fine.len(0), fine.len(1),
		                       nstencil, ibc);
	}


	void mpi_fortran_interp(const inter::ProlongOp & P,
	                        const core::GridFunc & coarse,
	                        const core::GridFunc & residual,
	                        core::GridFunc & fine)
	{
		using namespace boxmg::bmg2d;

		int nstencil, kf, kc, nog;

		inter::ProlongOp & Pd = const_cast<inter::ProlongOp&>(P);
		core::GridFunc & coarsed = const_cast<core::GridFunc&>(coarse);
		core::GridFunc & res = const_cast<core::GridFunc&>(residual);
		inter::mpi::ProlongOp & mpi_Pd = dynamic_cast<inter::mpi::ProlongOp&>(Pd);
		core::mpi::GridTopo & topo = mpi_Pd.grid();
		MsgCtx *ctx = (MsgCtx*) mpi_Pd.halo_ctx;

		if (Pd.stencil().five_pt()) {
			nstencil = 3;
		} else {
			nstencil = 5;
		}

		kc = topo.level() + 1;
		nog = topo.nlevel();
		kf = kc + 1;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_BMG2_SymStd_interp_add(kc, kf, nog,
		                           fine.data(), coarsed.data(), res.data(),
		                           mpi_Pd.fine_op->data(), nstencil,
		                           mpi_Pd.data(),
		                           coarsed.len(0), coarsed.len(1),
		                           fine.len(0), fine.len(1),
		                           topo.is(0), topo.is(1),
		                           ctx->msg_geom.data(), ctx->msg_geom.size(),
		                           ctx->pMSG.data(), ctx->msg_buffer.data(),
		                           ctx->msg_buffer.size(), fcomm);
	}
}


}}}
