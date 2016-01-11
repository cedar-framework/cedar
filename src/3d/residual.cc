#include <boxmg/types.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/mpi/halo.h>

#include <boxmg/3d/residual.h>

extern "C" {
	using namespace boxmg;
	void BMG3_SymStd_residual(int kg, int NOG, int ifd,
	                          real_t *q, real_t *qf, real_t *so, real_t *RES,
	                          len_t ii, len_t jj, len_t kk,
	                          int NStncl);
	void MPI_BMG3_SymStd_residual(int kg, int nog, int ifd,
	                              real_t *q, real_t *qf, real_t *so, real_t *res,
	                              len_t ii, len_t jj, len_t kk,
	                              int NStncl, len_t *iwork, len_t nmsgi,
	                              int *pMSG, real_t *msg_buffer, len_t nmsgr,
	                              int mpicomm);
}

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	void residual(const stencil_op & A, const grid_func & x,
	              const grid_func & b, grid_func & r)
	{
		int kg, nog, ifd, nstencil;

		auto &Ad = const_cast<stencil_op&>(A);
		auto &xd = const_cast<grid_func&>(x);
		auto &bd = const_cast<grid_func&>(b);

		if (Ad.stencil().five_pt()) {
			ifd = 1;
			nstencil = 4;
			nog = 1;
			kg = 1;
		} else {
			ifd = 0;
			nstencil = 14;
			nog = 1;
			kg = 1;
		}

		BMG3_SymStd_residual(kg, nog, ifd, xd.data(), bd.data(), Ad.data(), r.data(),
		                     r.len(0), r.len(1), r.len(2), nstencil);
	}


	void mpi_residual_fortran(const mpi::stencil_op & A, const mpi::grid_func & x,
	                          const mpi::grid_func & b, mpi::grid_func &r)
	{
		int k, kf, nog, ifd, nstencil;
		auto & Ad = const_cast<mpi::stencil_op &>(A);
		auto & xd = const_cast<mpi::grid_func&>(x);
		auto & bd = const_cast<mpi::grid_func&>(b);
		grid_topo & topo = Ad.grid();
		MsgCtx *ctx = (MsgCtx*) Ad.halo_ctx;

		if (Ad.stencil().five_pt()) {
			ifd = 1;
			nstencil = 4;
		} else {
			ifd = 0;
			nstencil = 14;
		}

		nog = kf = topo.nlevel();
		k = topo.level()+1;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		MPI_BMG3_SymStd_residual(k, nog, ifd,
		                         xd.data(), bd.data(), Ad.data(), r.data(),
		                         r.len(0), r.len(1), r.len(2),
		                         nstencil, ctx->msg_geom.data(), ctx->msg_geom.size(),
		                         ctx->pMSG.data(), ctx->msg_buffer.data(),
		                         ctx->msg_buffer.size(), fcomm);
	}
}

}}}
