#include <boxmg/3d/inter/setup_interp.h>
#include <boxmg/2d/ftn/BMG_parameters_c.h>
#include <boxmg/3d/mpi/halo.h>

extern "C" {
	using namespace boxmg;
	void BMG3_SymStd_SETUP_interp_OI(int kgf, int kgc, real_t *so, real_t *soc,
	                                 real_t *ci, len_t iif, len_t jjf, len_t kkf,
	                                 len_t iic, len_t jjc, len_t kkc,
	                                 int nog, int ifd, int nstncl, int irelax, int jpn,
	                                 real_t *yo);
	void MPI_BMG3_SymStd_SETUP_interp_OI(int kgf, int kgc, real_t *so, real_t *soc,
	                                     real_t *ci, len_t iif, len_t jjf, len_t kkf,
	                                     len_t iic, len_t jjc, len_t kkc,
	                                     int nog, int ifd, int nstencil, int irelax, real_t *yo,
	                                     int nogm, len_t *IGRD, len_t *iwork, len_t NMSGi,
	                                     int *pMSG, real_t *buffer, len_t nmsgr, int myproc,
	                                     int mpicomm);

}


namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	using namespace boxmg::bmg3;

	void setup_interp(int kf, int kc, int nog,
	                  const stencil_op & fop,
	                  const stencil_op & cop,
	                  inter::prolong_op &P)
	{
		int nstencil;
		int irelax = BMG_RELAX_SYM;
		int jpn = BMG_BCs_definite;
		int ifd;

		const grid_stencil &fsten = fop.stencil();
		const grid_stencil &csten = cop.stencil();
		stencil_op &fopd = const_cast<stencil_op&>(fop);
		stencil_op &copd = const_cast<stencil_op&>(cop);

		P.fine_op = &fopd;

		if (fsten.five_pt()) {
			ifd = 1;
			nstencil = 4;
		} else {
			ifd = 0;
			nstencil = 14;
		}

		// TODO: preallocate this?
		array<len_t, real_t, 4> yo(fsten.len(0), fsten.len(1), 2, 14);

		BMG3_SymStd_SETUP_interp_OI(kf, kc,
		                            fopd.data(), copd.data(),
		                            P.data(),
		                            fsten.len(0), fsten.len(1), fsten.len(2),
		                            csten.len(0), csten.len(1), csten.len(2),
		                            nog, ifd, nstencil, irelax, jpn,
		                            yo.data());
	}


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
			nstencil = 4;
		} else {
			ifd = 0;
			nstencil = 14;
		}

		int rank;
		MPI_Comm_rank(topo.comm, &rank);
		rank++;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		array<len_t, real_t, 4> yo(fsten.len(0), fsten.len(1), 2, 14);
		MPI_BMG3_SymStd_SETUP_interp_OI(kf, kc, fopd.data(), copd.data(),
		                                P.data(), fsten.len(0), fsten.len(1), fsten.len(2),
		                                csten.len(0), csten.len(1), csten.len(2),
		                                nog, ifd, nstencil, BMG_RELAX_SYM,
		                                yo.data(), nog, topo.IGRD(),
		                                ctx->msg_geom.data(), ctx->msg_geom.size(),
		                                ctx->pMSG.data(), ctx->msg_buffer.data(), ctx->msg_buffer.size(),
		                                rank, fcomm);

	}
}

}}}
