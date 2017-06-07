#ifndef CEDAR_3D_HALO_H
#define CEDAR_3D_HALO_H

#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"

#include <cedar/kernel_params.h>
#include <cedar/array.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/grid_func.h>

namespace cedar { namespace cdr3 { namespace kernel {
namespace impls
{
	namespace mpi = cedar::cdr3::mpi;
	struct MsgCtx
	{
		MsgCtx(grid_topo & topo);
		array<int,int,2> pMSG; // these should be int, len_t
		array<int,int,2> pLS;
		array<int,int,2> pMSGSO;
		std::vector<len_t> msg_geom;
		array<int,int,3> proc_grid;
		std::vector<int> proc_coord;
		std::vector<len_t> dimxfine;
		std::vector<len_t> dimyfine;
		std::vector<len_t> dimzfine;
		array<int,len_t,2> dimx;
		array<int,len_t,2> dimy;
		array<int,len_t,2> dimz;
		std::vector<real_t> msg_buffer;
		int pSI_MSG;
		int p_NLx_kg, p_NLy_kg, p_NLz_kg;
		/* std::vector<len_t> iworkmsg; */
		/* int *iworkmsg[nmsgi]; */
		/* int nmsgi; */
		/* int pmsg[nbmg_pmsg,nog]; */
		/* int msgbuffer[nmsgr]; */
		/* int nmsgr; */
		len_t nlocal(int kg, int dim, int ijkrank) {
			len_t local_arr_ptr = pMSG(ipL_MSG_LocalArraySize,kg) - 1;  // 1 vs 0 based indexing
			len_t idx = local_arr_ptr;
			ijkrank--; // 1 vs 0 based indexing
			idx += ijkrank*3 + dim;

			return msg_geom[idx];
		}

		len_t cg_nlocal(int dim, int ijkrank) {
			return nlocal(0, dim, ijkrank);
		}
	};

	void setup_msg(const kernel_params & params, grid_topo &topo, void **msg_ctx);
	void msg_exchange(const kernel_params & params, mpi::grid_func & f);
	template<class sten>
	void msg_stencil_exchange(const kernel_params & params, mpi::stencil_op<sten> & so)
	{
		MsgCtx *ctx = (MsgCtx*) sop.halo_ctx;
		grid_topo &topo = sop.grid();
		int nstencil = stencil_ndirs<sten>::value;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		BMG3_SymStd_SETUP_fine_sopcil(topo.level()+1, sop.data(),
		                               sop.len(0), sop.len(1), sop.len(2),
		                               nstencil,
		                               ctx->msg_geom.data(), ctx->msg_geom.size(),
		                               ctx->pMSGSO.data(), ctx->msg_buffer.data(),
		                               ctx->msg_buffer.size(), fcomm);
	}
}
}}}

#endif
