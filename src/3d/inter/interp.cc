#include <boxmg/2d/ftn/BMG_parameters_c.h>
#include <boxmg/3d/mpi/halo.h>

#include <boxmg/3d/inter/interp.h>

extern "C" {
	using namespace boxmg;
	void BMG3_SymStd_interp_add(real_t *q, real_t *qc,
	                            real_t *so, real_t *res, real_t *ci,
	                            len_t iic, len_t jjc, len_t kkc,
	                            len_t iif, len_t jjf, len_t kkf,
	                            int NStncl, int jpn);
	void MPI_BMG3_SymStd_interp_add(int kcg, int kfg, int nog,
	                                real_t *q, real_t *qc, real_t *res, real_t *so,
	                                int nstncl, real_t *ci,
	                                len_t iic, len_t jjc, len_t kkc,
	                                len_t iif, len_t jjf, len_t kkf,
	                                len_t igs, len_t jgs, len_t kgs,
	                                len_t *iwork, len_t nmsgi, int *pmsg,
	                                real_t *buffer, len_t nmsgr, int mpicomm);
}

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void fortran_interp(const inter::prolong_op & P,
	                    const grid_func & coarse,
	                    const grid_func & residual,
	                    grid_func & fine)
	{
		using namespace boxmg::bmg3;

		int nstencil, ibc;

		inter::prolong_op & Pd = const_cast<inter::prolong_op&>(P);
		grid_func & coarsed = const_cast<grid_func&>(coarse);
		grid_func & res = const_cast<grid_func&>(residual);

		if (Pd.stencil().five_pt()) {
			nstencil = 4;
		} else {
			nstencil = 14;
		}

		ibc = BMG_BCs_definite;

		BMG3_SymStd_interp_add(fine.data(), coarsed.data(), Pd.fine_op->data(), res.data(), Pd.data(),
		                       coarsed.len(0), coarsed.len(1), coarsed.len(2),
		                       fine.len(0), fine.len(1), fine.len(2),
		                       nstencil, ibc);
	}


	void mpi_fortran_interp(const inter::mpi::prolong_op & P,
	                        const mpi::grid_func & coarse,
	                        const mpi::grid_func & residual,
	                        mpi::grid_func & fine)
	{
		using namespace boxmg::bmg3;

		int nstencil, kf, kc, nog;

		inter::mpi::prolong_op & Pd = const_cast<inter::mpi::prolong_op&>(P);
		mpi::grid_func & coarsed = const_cast<mpi::grid_func&>(coarse);
		mpi::grid_func & res = const_cast<mpi::grid_func&>(residual);
		grid_topo & topo = Pd.grid();
		MsgCtx *ctx = (MsgCtx*) Pd.halo_ctx;

		if (Pd.stencil().five_pt()) {
			nstencil = 4;
		} else {
			nstencil = 14;
		}

		kc = topo.level() + 1;
		nog = topo.nlevel();
		kf = kc + 1;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_BMG3_SymStd_interp_add(kc, kf, nog,
		                           fine.data(), coarsed.data(), res.data(),
		                           Pd.fine_op->data(), nstencil, Pd.data(),
		                           coarsed.len(0), coarsed.len(1), coarsed.len(2),
		                           fine.len(0), fine.len(1), fine.len(2),
		                           topo.is(0), topo.is(1), topo.is(2),
		                           ctx->msg_geom.data(), ctx->msg_geom.size(),
		                           ctx->pMSG.data(), ctx->msg_buffer.data(),
		                           ctx->msg_buffer.size(), fcomm);
	}
}

}}}
