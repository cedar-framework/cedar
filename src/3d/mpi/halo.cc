#include <cstdlib>
#include <iostream>

#include <cedar/types.h>
#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>
#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>

#include <cedar/3d/mpi/halo.h>
#include <cedar/3d/mpi/msg_exchanger.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_MSG(int *pMSG, int *pMSGSO, len_t *imsg_geom,
	                           len_t nmsgi, int *pSI_MSG, int IBC, len_t *IGRD,
	                           int nog, int nogm, int nproc, int myproc,
	                           len_t *dimx, len_t *dimy, len_t *dimz,
	                           len_t *dimxfine, len_t *dimyfine, len_t *dimzfine,
	                           int *proc_grid, int nproci, int nprocj, int nprock,
	                           int mpicomm);
	void BMG3_SymStd_SETUP_fine_stencil(int kf, real_t *so,
	                                    len_t nlx, len_t nly, len_t nlz,
	                                    int nstencil,
	                                    len_t *iwork, len_t nmsgi, int *pMSGSO,
	                                    real_t *buffer, len_t nmsgr,
	                                    int mpicomm);
	void BMG3_SymStd_UTILS_update_ghosts(int k, real_t *x,
	                                     len_t nx, len_t ny, len_t nz,
	                                     len_t *iwork, len_t nmsgi, int *pMSG,
	                                     real_t *buffer, len_t nmsgr, int nog,
	                                     int mpicomm);
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr3 { namespace kernel {
namespace impls
{
	std::unique_ptr<halo_exchanger> setup_msg(const kernel_params & params, grid_topo & topo)
	{
		auto halof = std::make_unique<msg_exchanger>(topo);
		auto & ctx = halof->context();
		int rank;
		int ibc;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		MPI_Comm_rank(topo.comm, &rank);
		rank++; // Fortran likes to be difficult...

		BMG_get_bc(params.per_mask(), &ibc);

		BMG3_SymStd_SETUP_MSG(ctx.pMSG.data(), ctx.pMSGSO.data(),
		                      ctx.msg_geom.data(), ctx.msg_geom.size(),
		                      &ctx.pSI_MSG, ibc, topo.IGRD(), topo.nlevel(),
		                      topo.nlevel(), topo.nproc(), rank,
		                      ctx.dimx.data(), ctx.dimy.data(), ctx.dimz.data(),
		                      ctx.dimxfine.data(), ctx.dimyfine.data(), ctx.dimzfine.data(),
		                      ctx.proc_grid.data(), topo.nproc(0), topo.nproc(1), topo.nproc(2),
		                      fcomm);

		halof->init();
		return halof;
	}

	void msg_exchange(const kernel_params & params, halo_exchanger *halof, mpi::grid_func & f)
	{
		MsgCtx *ctx = (MsgCtx*) halof->context_ptr();
		grid_topo &topo = f.grid();

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		BMG3_SymStd_UTILS_update_ghosts(topo.level()+1, f.data(),
		                                f.len(0), f.len(1), f.len(2),
		                                ctx->msg_geom.data(), ctx->msg_geom.size(),
		                                ctx->pMSG.data(), ctx->msg_buffer.data(),
		                                ctx->msg_buffer.size(), topo.nlevel(), fcomm);
	}
}
}}}
