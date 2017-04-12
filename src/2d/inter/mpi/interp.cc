#include "cedar/2d/ftn/BMG_parameters_c.h"
#include "cedar/2d/inter/mpi/prolong_op.h"
#include "cedar/2d/mpi/halo.h"

#include "cedar/2d/inter/interp.h"

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_interp_add(len_t kc, len_t kf, len_t nog,
	                                real_t *Q, real_t *QC, real_t *RES, real_t *SO,
	                                int nstncl, real_t *CI, len_t iic, len_t jjc,
	                                len_t iif, len_t jjf, len_t iGs, len_t jGs,
	                                len_t *iWorkMSG, len_t NMSGi, int *pMSG,
	                                real_t *msg_buffer, len_t nmsgr, int mpicomm);
}

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	void mpi_fortran_interp(const kernel_params & params,
	                        const inter::mpi::prolong_op & P,
	                        const mpi::grid_func & coarse,
	                        const mpi::grid_func & residual,
	                        mpi::grid_func & fine)
	{
		using namespace cedar::cdr2;

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
	}
}


}}}
