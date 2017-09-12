#ifndef CEDAR_2D_KERNEL_RESIDUAL_H
#define CEDAR_2D_KERNEL_RESIDUAL_H

#include <cedar/kernel_params.h>
#include <cedar/2d/mpi/msg_exchanger.h>
#include <cedar/2d/grid_func.h>
#include <cedar/2d/stencil_op.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/grid_func.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_residual(int k, int kf, int nog,
	                              real_t *SO, real_t *QF, real_t *Q, real_t *RES,
	                              len_t ii, len_t jj, int ifd, int nstencil,
	                              int irelax, int irelax_sym,
	                              len_t *iWorkMSG, len_t NMSGi, int *pMSG,
	                              real_t *msg_buffer, len_t nmsgr, int mpicomm);
}


namespace cedar { namespace cdr2 { namespace kernel {
namespace impls
{
	template<class sten>
	void residual_fortran(const kernel_params & params, const stencil_op<sten> & A, const grid_func & x,
						  const grid_func & b, grid_func &y);

	namespace mpi = cedar::cdr2::mpi;
	template <class sten>
		void mpi_residual_fortran(const kernel_params & params,
		                          mpi::msg_exchanger *halof,
		                          const mpi::stencil_op<sten> & A,
		                          const mpi::grid_func & x,
		                          const mpi::grid_func & b, mpi::grid_func &r)
	{
		int k, kf, nog, ifd, nstencil;
		auto & Ad = const_cast<mpi::stencil_op<sten> &>(A);
		auto & xd = const_cast<mpi::grid_func&>(x);
		auto & bd = const_cast<mpi::grid_func&>(b);
		grid_topo & topo = Ad.grid();
		MsgCtx *ctx = (MsgCtx*) halof->context_ptr();

		nstencil = stencil_ndirs<sten>::value;
		if (std::is_same<five_pt, sten>::value)
			ifd = 1;
		else
			ifd = 0;

		int irelax = 0;
		int irelax_sym = 0;

		nog = kf = topo.nlevel();
		k = topo.level()+1;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_BMG2_SymStd_residual(k, kf, nog,
		                         Ad.data(), bd.data(), xd.data(), r.data(),
		                         r.len(0), r.len(1), ifd, nstencil,
		                         irelax, irelax_sym,
		                         ctx->msg_geom.data(), ctx->msg_geom.size(),
		                         ctx->pMSG.data(), ctx->msg_buffer.data(),
		                         ctx->msg_buffer.size(), fcomm);
		halof->exchange_func(k, r.data());
	}
}

}}}

#endif
