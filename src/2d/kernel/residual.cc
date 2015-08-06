#include "core/mpi/stencil_op.h"
#include "halo.h"

#include "residual.h"

	extern "C" {
		using namespace boxmg;
		void BMG2_SymStd_residual(int*,real_t*,real_t*,real_t*,real_t*,len_t*,len_t*,
								  int*,int*,int*,int*,int*,int*,int*);
		void MPI_BMG2_SymStd_residual(int k, int kf, int nog,
		                              real_t *SO, real_t *QF, real_t *Q, real_t *RES,
		                              len_t ii, len_t jj, int ifd, int nstencil,
		                              int irelax, int irelax_sym,
		                              len_t *iWorkMSG, len_t NMSGi, int *pMSG,
		                              real_t *msg_buffer, len_t nmsgr, int mpicomm);
	}

namespace boxmg { namespace bmg2d { namespace kernel {


namespace impls
{
	using namespace boxmg::bmg2d;

	void residual(const core::StencilOp & A, const core::GridFunc & x,
				  const core::GridFunc & b, core::GridFunc &r)
	{
		using namespace boxmg::bmg2d::core;

		const GridStencil &so = A.stencil();

		for (auto j: r.range(1)) {
			for (auto i: r.range(0)) {
				r(i,j) = (b(i,j) +
				          (so(i,j,Dir::W)  * x(i-1, j  ) +
				           so(i,j,Dir::E)  * x(i+1, j  ) +
				           so(i,j,Dir::S)  * x(i  , j-1) +
				           so(i,j,Dir::N)  * x(i  , j+1) +
						   so(i,j,Dir::SW) * x(i-1, j-1) +
				           so(i,j,Dir::SE) * x(i+1, j-1) +
				           so(i,j,Dir::NW) * x(i-1, j+1) +
				           so(i,j,Dir::NE) * x(i+1, j+1) -
						   so(i,j,Dir::C)  * x(i  , j)));
			}
		}
	}

	void residual_fortran(const core::StencilOp &A, const core::GridFunc &x,
						  const core::GridFunc &b, core::GridFunc &r)
	{
		int k = 0;
		int kf = 0;
		int ifd = 0;
		int nstncl = 5;
		int ibc = 0;
		int irelax = 0;
		int irelax_sym = 0;
		int updown = 0;
		len_t ii = r.len(0);
		len_t jj = r.len(1);

		core::StencilOp &Ad = const_cast<core::StencilOp&>(A);
		core::GridFunc &xd = const_cast<core::GridFunc&>(x);
		core::GridFunc &bd = const_cast<core::GridFunc&>(b);
		using namespace boxmg::bmg2d::core;

		BMG2_SymStd_residual(&k, Ad.data(), bd.data(), xd.data(), r.data(), &ii, &jj,
							 &kf, &ifd, &nstncl, &ibc, &irelax, &irelax_sym, &updown);

	}


	void mpi_residual_fortran(const core::StencilOp &A, const core::GridFunc &x,
	                          const core::GridFunc &b, core::GridFunc &r)
	{
		int k, kf, nog, ifd, nstencil;
		auto & Ad = const_cast<core::StencilOp &>(A);
		auto & xd = const_cast<core::GridFunc&>(x);
		auto & bd = const_cast<core::GridFunc&>(b);
		auto & mpi_Ad = dynamic_cast<core::mpi::StencilOp&>(Ad);
		core::mpi::GridTopo & topo = mpi_Ad.grid();
		MsgCtx *ctx = (MsgCtx*) mpi_Ad.halo_ctx;

		if (Ad.stencil().five_pt()) {
			ifd = 1;
			nstencil = 3;
		} else {
			ifd = 0;
			nstencil = 5;
		}

		int irelax = 0;
		int irelax_sym = 0;

		nog = kf = topo.nlevel();
		k = topo.level()+1;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_BMG2_SymStd_residual(k, kf, nog,
		                         mpi_Ad.data(), bd.data(), xd.data(), r.data(),
		                         r.len(0), r.len(1), ifd, nstencil,
		                         irelax, irelax_sym,
		                         ctx->msg_geom.data(), ctx->msg_geom.size(),
		                         ctx->pMSG.data(), ctx->msg_buffer.data(),
		                         ctx->msg_buffer.size(), fcomm);
	}
}


}}}
