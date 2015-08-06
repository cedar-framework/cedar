#ifndef BOXMG_2D_KERNEL_HALO_H
#define BOXMG_2D_KERNEL_HALO_H

#include "core/mpi/grid_topo.h"
#include "core/mpi/stencil_op.h"
#include "core/mpi/grid_func.h"
#include "core/array.h"

namespace boxmg { namespace bmg2d { namespace kernel {
namespace impls
{
	struct MsgCtx
	{
		MsgCtx(core::mpi::GridTopo & topo);
		core::Array<int,int> pMSG; // these should be int, len_t
		core::Array<int,int> pMSGSO;
		std::vector<len_t> msg_geom;
		core::Array<int,int> proc_grid;
		std::vector<int> proc_coord;
		std::vector<int> dimxfine;
		std::vector<int> dimyfine;
		core::Array<int,int> dimx;
		core::Array<int,int> dimy;
		std::vector<real_t> msg_buffer;
		int pSI_MSG;
		int p_NLx_kg, p_NLy_kg;
		/* std::vector<len_t> iworkmsg; */
		/* int *iworkmsg[nmsgi]; */
		/* int nmsgi; */
		/* int pmsg[nbmg_pmsg,nog]; */
		/* int msgbuffer[nmsgr]; */
		/* int nmsgr; */
	};

	void setup_msg(core::mpi::GridTopo &topo, void **msg_ctx);
	void msg_exchange(core::mpi::GridFunc & f);
	void msg_stencil_exchange(core::mpi::StencilOp & so);
}
}}}

#endif
