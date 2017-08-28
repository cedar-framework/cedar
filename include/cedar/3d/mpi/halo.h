#ifndef CEDAR_3D_HALO_H
#define CEDAR_3D_HALO_H

#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"

#include <cedar/kernel_params.h>
#include <cedar/halo_exchanger.h>
#include <cedar/array.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/grid_func.h>
#include <cedar/3d/mpi/msg_exchanger.h>


extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_fine_stencil(int kf, real_t *so,
	                                    len_t nlx, len_t nly, len_t nlz,
	                                    int nstencil,
	                                    len_t *iwork, len_t nmsgi, int *pMSGSO,
	                                    real_t *buffer, len_t nmsgr,
	                                    int mpicomm);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr3 { namespace kernel {
namespace impls
{
	namespace mpi = cedar::cdr3::mpi;

	std::unique_ptr<halo_exchanger> setup_msg(const kernel_params & params, grid_topo &topo);
	void msg_exchange(const kernel_params & params, halo_exchanger *halof, mpi::grid_func & f);
	template<class sten>
		void msg_stencil_exchange(const kernel_params & params, halo_exchanger *halof, mpi::stencil_op<sten> & sop)
	{
		MsgCtx *ctx = (MsgCtx*) halof->context_ptr();
		grid_topo &topo = sop.grid();
		int nstencil = stencil_ndirs<sten>::value;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		BMG3_SymStd_SETUP_fine_stencil(topo.level()+1, sop.data(),
		                               sop.len(0), sop.len(1), sop.len(2),
		                               nstencil,
		                               ctx->msg_geom.data(), ctx->msg_geom.size(),
		                               ctx->pMSGSO.data(), ctx->msg_buffer.data(),
		                               ctx->msg_buffer.size(), fcomm);
	}
}
}}}

#endif
