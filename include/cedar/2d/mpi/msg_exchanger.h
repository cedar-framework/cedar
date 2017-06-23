#ifndef CEDAR_MSG_EXCHANGER_H
#define CEDAR_MSG_EXCHANGER_H

#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"

#include <cedar/halo_exchanger.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/grid_func.h>

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
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
	MPI_Comm comm;
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

class msg_exchanger : public halo_exchanger
{
public:
	msg_exchanger(grid_topo & topo);
	MsgCtx & context() { return ctx; }
	void *context_ptr() { return &ctx; }
	virtual void exchange_func(int k, real_t *gf) override;
	virtual void exchange_sten(int k, real_t *so) override;

private:
	MsgCtx ctx;
	/**
	   Local array extents by grid number needed for MSG ghost
	   exchange calls.  The first dimension is the grid number,
	   the second is the dimension.
	**/
	array<len_t, len_t, 2> dims;
};


}}}}

#endif
