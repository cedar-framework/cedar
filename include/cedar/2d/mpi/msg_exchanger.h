#ifndef CEDAR_2D_MSG_EXCHANGER_H
#define CEDAR_2D_MSG_EXCHANGER_H

#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>
#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>

#include <cedar/services/halo_exchange.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/2d/mpi/types.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/grid_func.h>

namespace cedar { namespace cdr2 {

namespace impls
{
struct MsgCtx
{
	MsgCtx(grid_topo & topo);
	aarray<int,int,2> pMSG; // these should be int, len_t
	aarray<int,int,2> pLS;
	aarray<int,int,2> pMSGSO;
	std::vector<len_t> msg_geom;
	aarray<int,int,2> proc_grid;
	std::vector<int> proc_coord;
	std::vector<len_t> dimxfine;
	std::vector<len_t> dimyfine;
	aarray<int,len_t,2> dimx;
	aarray<int,len_t,2> dimy;
	std::vector<real_t> msg_buffer;
	int pSI_MSG;
	int p_NLx_kg, p_NLy_kg;
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
}

namespace mpi {

class msg_exchanger : public services::halo_exchange<stypes>
{
	using MsgCtx = impls::MsgCtx;
public:
	using parent = services::halo_exchange<stypes>;
	MsgCtx & context() { return *ctx; }
	void *context_ptr() { return ctx.get(); }
	void setup(std::vector<topo_ptr> topos) override;
	void run(stencil_op<five_pt> & so) override;
	void run(stencil_op<nine_pt> & so) override;
	void run(grid_func & gf, unsigned short dmask=7) override;
	void exchange_func(int k, real_t *gf) override;
	void exchange_sten(int k, real_t *so) override;
	aarray<int, len_t, 2> & leveldims(int k) override {
		if (k == 0)
			return ctx->dimx;
		else
			return ctx->dimy;
	}
	len_t * datadist(int k, int grid) override {
		if (k == 0)
			return &ctx->msg_geom.data()[ctx->pLS(ipL_LS_XDataDist,grid-1)-1];
		else
			return &ctx->msg_geom.data()[ctx->pLS(ipL_LS_YDataDist,grid-1)-1];
	}
	std::vector<real_t> & linebuf() override {
		return ctx->msg_buffer;
	}

private:
	std::unique_ptr<MsgCtx> ctx;
	/**
	   Local array extents by grid number needed for MSG ghost
	   exchange calls.  The first dimension is the grid number,
	   the second is the dimension.
	**/
	array<len_t, 2> dims;
	array<len_t, 1> coord;
};


}}}

#endif
