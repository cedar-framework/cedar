#ifndef BOXMG_3D_HALO_H
#define BOXMG_3D_HALO_H

#include "boxmg/2d/ftn/mpi/BMG_workspace_c.h"

#include <boxmg/array.h>
#include <boxmg/mpi/grid_topo.h>
#include <boxmg/3d/mpi/stencil_op.h>
#include <boxmg/3d/mpi/grid_func.h>

namespace boxmg { namespace bmg3 { namespace kernel {
namespace impls
{
	namespace mpi = boxmg::bmg3::mpi;
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

	void setup_msg(grid_topo &topo, void **msg_ctx);
	void msg_exchange(mpi::grid_func & f);
	void msg_stencil_exchange(mpi::stencil_op & so);
}
}}}

#endif
