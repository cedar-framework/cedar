#include "boxmg/2d/ftn/BMG_parameters_c.h"
#include "boxmg/2d/inter/mpi/prolong_op.h"
#include "boxmg/2d/mpi/halo.h"
#include <boxmg/util/timer.h>

#include "boxmg/2d/inter/interp.h"

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
	void fortran_interp(const inter::prolong_op & P,
	                    const grid_func & coarse,
	                    const grid_func & residual,
	                    grid_func & fine)
	{
		using namespace boxmg::bmg2d;

		int nstencil, ibc;

		inter::prolong_op & Pd = const_cast<inter::prolong_op&>(P);
		grid_func & coarsed = const_cast<grid_func&>(coarse);
		grid_func & res = const_cast<grid_func&>(residual);

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


	void mpi_fortran_interp(const inter::mpi::prolong_op & P,
	                        const mpi::grid_func & coarse,
	                        const mpi::grid_func & residual,
	                        mpi::grid_func & fine)
	{
		using namespace boxmg::bmg2d;

		int nstencil, kf, kc, nog;

		inter::mpi::prolong_op & Pd = const_cast<inter::mpi::prolong_op&>(P);
		mpi::grid_func & coarsed = const_cast<mpi::grid_func&>(coarse);
		mpi::grid_func & res = const_cast<mpi::grid_func&>(residual);
		grid_topo & topo = Pd.grid();
		grid_topo & topof = Pd.fine_op->grid();
		MsgCtx *ctx = (MsgCtx*) Pd.halo_ctx;

		if (Pd.stencil().five_pt()) {
			nstencil = 3;
		} else {
			nstencil = 5;
		}

		kc = topo.level() + 1;
		nog = topo.nlevel();
		kf = kc + 1;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		sync_timer itimer(topo.comm, "Interpolation");
		itimer.begin();
		MPI_BMG2_SymStd_interp_add(kc, kf, nog,
		                           fine.data(), coarsed.data(), res.data(),
		                           Pd.fine_op->data(), nstencil,
		                           Pd.data(),
		                           coarsed.len(0), coarsed.len(1),
		                           fine.len(0), fine.len(1),
		                           topof.is(0), topof.is(1),
		                           ctx->msg_geom.data(), ctx->msg_geom.size(),
		                           ctx->pMSG.data(), ctx->msg_buffer.data(),
		                           ctx->msg_buffer.size(), fcomm);
		itimer.end();
	}
}


}}}
