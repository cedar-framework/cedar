#include <cedar/2d/ftn/BMG_parameters_c.h>
#include <cedar/3d/mpi/halo.h>

#include <cedar/3d/inter/galerkin_prod.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_SETUP_ITLI07_ex(int kgf, int kgc, real_t *so, real_t *soc, real_t *ci,
	                                     len_t iif, len_t jjf, len_t kkf,
	                                     len_t iic, len_t jjc, len_t kkc,
	                                     len_t iGs, len_t jGs, len_t kGs,
	                                     int nog, int nogm, len_t *IGRD,
	                                     len_t *iwork, len_t NMSGi,
	                                     int *pMSGSO, real_t *buffer, len_t NMSGr,
	                                     int nproc, int *proc_grid, int myproci,
	                                     int myprocj, int myprock, int nproci, int nprocj,
	                                     int nprock, len_t *dimx, len_t *dimy, len_t *dimz, int mpicomm);
	void MPI_BMG3_SymStd_SETUP_ITLI27_ex(int kgf, int kgc, real_t *so, real_t *soc, real_t *ci,
	                                     len_t iif, len_t jjf, len_t kkf,
	                                     len_t iic, len_t jjc, len_t kkc,
	                                     len_t iGs, len_t jGs, len_t kGs,
	                                     int nog, int nogm, len_t *IGRD,
	                                     len_t *iwork, len_t NMSGi, int *pSI_MSG,
	                                     int *pMSGSO, real_t *buffer, len_t NMSGr,
	                                     int nproc, int *proc_grid, int myproci,
	                                     int myprocj, int myprock, int nproci, int nprocj,
	                                     int nprock, len_t *dimx, len_t *dimy, len_t *dimz, int mpicomm);
}

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr3;

	void mpi_galerkin_prod(const kernel_params & params,
	                       int kf, int kc, int nog,
	                       const inter::mpi::prolong_op & P,
	                       const mpi::stencil_op & fop,
	                       mpi::stencil_op & copd)
	{
		inter::mpi::prolong_op & Pd = const_cast<inter::mpi::prolong_op&>(P);
		mpi::stencil_op & fopd = const_cast<mpi::stencil_op&>(fop);
		grid_topo & topo = fopd.grid();
		grid_stencil & fsten = fopd.stencil();
		grid_stencil & csten = copd.stencil();
		MsgCtx *ctx = (MsgCtx*) fopd.halo_ctx;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		if (fsten.five_pt()) {
			MPI_BMG3_SymStd_SETUP_ITLI07_ex(kf, kc, fopd.data(), copd.data(), Pd.data(),
			                                fsten.len(0), fsten.len(1), fsten.len(2),
			                                csten.len(0), csten.len(1), csten.len(2),
			                                topo.is(0), topo.is(1), topo.is(2),
			                                nog, nog, topo.IGRD(),
			                                ctx->msg_geom.data(), ctx->msg_geom.size(),
			                                ctx->pMSGSO.data(),
			                                ctx->msg_buffer.data(), ctx->msg_buffer.size(),
			                                topo.nproc(), ctx->proc_grid.data(),
			                                topo.coord(0), topo.coord(1), topo.coord(2),
			                                topo.nproc(0), topo.nproc(1), topo.nproc(2),
			                                ctx->dimx.data(), ctx->dimy.data(), ctx->dimz.data(),
			                                fcomm);
		} else {
			MPI_BMG3_SymStd_SETUP_ITLI27_ex(kf, kc, fopd.data(), copd.data(), Pd.data(),
			                                fsten.len(0), fsten.len(1), fsten.len(2),
			                                csten.len(0), csten.len(1), csten.len(2),
			                                topo.is(0), topo.is(1), topo.is(2),
			                                nog, nog, topo.IGRD(),
			                                ctx->msg_geom.data(), ctx->msg_geom.size(),
			                                &ctx->pSI_MSG, ctx->pMSGSO.data(),
			                                ctx->msg_buffer.data(), ctx->msg_buffer.size(),
			                                topo.nproc(), ctx->proc_grid.data(),
			                                topo.coord(0), topo.coord(1), topo.coord(2),
			                                topo.nproc(0), topo.nproc(1), topo.nproc(2),
			                                ctx->dimx.data(), ctx->dimy.data(), ctx->dimz.data(),
			                                fcomm);
		}


	}
}

}}}
