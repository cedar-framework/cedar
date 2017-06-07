#ifndef CEDAR_3D_KERNEL_MATVEC_H
#define CEDAR_3D_KERNEL_MATVEC_H

#include <type_traits>

#include <cedar/kernel_params.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/grid_func.h>
#include <cedar/3d/mpi/halo.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_UTILS_matvec(int kg, real_t *so, real_t *qf,
	                                  real_t *q, len_t ii, len_t jj, len_t kk,
	                                  int nog, int ifd, int nstncl,
	                                  len_t *iwork, int *pMSG, real_t *msg_buffer,
	                                  int mpicomm);
}

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;
	template<class sten>
	void matvec(const kernel_params & params,
	            const mpi::stencil_op<sten> & so,
	            const mpi::grid_func & x,
	            mpi::grid_func & y)
	{
		using namespace cedar::cdr3;
		int kg, ifd, nstencil, nog;

		auto & sod = const_cast<mpi::stencil_op<sten>&>(so);
		mpi::grid_func & xd = const_cast<mpi::grid_func&>(x);
		grid_topo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) sod.halo_ctx;

		nog = topo.nlevel();
		kg = topo.level()+1;
		nstencil = stencil_ndirs<sten>::value;
		if (std::is_same<sten, seven_pt>::value)
			ifd = 1;
		else
			ifd = 0;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		MPI_BMG3_SymStd_UTILS_matvec(kg, sod.data(), b.data(), xd.data(),
		                             so.len(0), so.len(1), so.len(2),
		                             nog, ifd, nstencil,
		                             ctx->msg_geom.data(), ctx->pMSG.data(),
		                             ctx->msg_buffer.data(), fcomm);
	}
}

}}}
#endif
