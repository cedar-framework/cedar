#ifndef BOXMG_2D_KERNEL_HALO_H
#define BOXMG_2D_KERNEL_HALO_H

#include "boxmg/2d/ftn/mpi/BMG_workspace_c.h"

#include <boxmg/array.h>
#include <boxmg/mpi/grid_topo.h>
#include <boxmg/2d/mpi/stencil_op.h>
#include <boxmg/2d/mpi/grid_func.h>

namespace boxmg { namespace bmg2d { namespace kernel {
namespace impls
{
	namespace mpi = boxmg::bmg2d::mpi;
	struct MsgCtx
	{
		MsgCtx(grid_topo & topo);
		array<int,int,2> pMSG; // these should be int, len_t
		array<int,int,2> pLS;
		array<int,int,2> pMSGSO;
		std::vector<len_t> msg_geom;
		array<int,int,2> proc_grid;
		std::vector<int> proc_coord;
		std::vector<len_t> dimxfine;
		std::vector<len_t> dimyfine;
		array<int,len_t,2> dimx;
		array<int,len_t,2> dimy;
		std::vector<real_t> msg_buffer;
		int pSI_MSG;
		int p_NLx_kg, p_NLy_kg;
		MPI_Comm xlinecomm;
		MPI_Comm ylinecomm;
		/* std::vector<len_t> iworkmsg; */
		/* int *iworkmsg[nmsgi]; */
		/* int nmsgi; */
		/* int pmsg[nbmg_pmsg,nog]; */
		/* int msgbuffer[nmsgr]; */
		/* int nmsgr; */
		len_t nlocal(int kg, int dim, int ijrank) {
			len_t local_arr_ptr = pMSG(ipL_MSG_LocalArraySize,kg) - 1;  // 1 vs 0 based indexing
			len_t idx = local_arr_ptr;
			ijrank--; // 1 vs 0 based indexing
			idx += ijrank*3 + dim;

			return msg_geom[idx];
		}

		len_t cg_nlocal(int dim, int ijrank) {
			return nlocal(0, dim, ijrank);
		}
	};

	void setup_msg(grid_topo &topo, void **msg_ctx);
	void msg_exchange(mpi::grid_func & f);
	void msg_stencil_exchange(mpi::stencil_op & so);
}
}}}

#endif
