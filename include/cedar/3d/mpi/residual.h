#ifndef CEDAR_3D_KERNEL_RESIDUAL_H
#define CEDAR_3D_KERNEL_RESIDUAL_H

#include <type_traits>
#include <cedar/kernel_params.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/grid_func.h>
#include <cedar/3d/mpi/halo.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_residual(int kg, int nog, int ifd,
	                              real_t *q, real_t *qf, real_t *so, real_t *res,
	                              len_t ii, len_t jj, len_t kk,
	                              int NStncl, len_t *iwork, len_t nmsgi,
	                              int *pMSG, real_t *msg_buffer, len_t nmsgr,
	                              int mpicomm);
}

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;
	template<class sten>
	void mpi_residual_fortran(const kernel_params & params, const mpi::stencil_op<sten> & A, const mpi::grid_func & x,
	                          const mpi::grid_func & b, mpi::grid_func &y)
	{
		int k, kf, nog, ifd, nstencil;
		auto & Ad = const_cast<mpi::stencil_op<sten> &>(A);
		auto & xd = const_cast<mpi::grid_func&>(x);
		auto & bd = const_cast<mpi::grid_func&>(b);
		grid_topo & topo = Ad.grid();
		MsgCtx *ctx = (MsgCtx*) Ad.halo_ctx;

		nstencil = stencil_ndirs<sten>::value;
		if (std::is_same<sten, seven_pt>::value)
			ifd = 1;
		else
			ifd = 0;

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
#endif
