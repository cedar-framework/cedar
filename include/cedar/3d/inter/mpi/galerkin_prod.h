#ifndef CEDAR_3D_INTER_MPI_GALERKIN_PROD_H
#define CEDAR_3D_INTER_MPI_GALERKIN_PROD_H

#include <type_traits>

#include <cedar/kernel_params.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/halo.h>
#include <cedar/3d/inter/mpi/prolong_op.h>
#include <cedar/2d/ftn/BMG_parameters_c.h>

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
	namespace mpi = cedar::cdr3::mpi;

	template<class sten>
	void mpi_galerkin_prod(const kernel_params & params,
	                       const inter::mpi::prolong_op & P,
	                       const mpi::stencil_op<sten> & fop,
	                       mpi::stencil_op<xxvii_pt> & cop)
	{
		int kf, kc, nog;
		inter::mpi::prolong_op & Pd = const_cast<inter::mpi::prolong_op&>(P);
		auto & fopd = const_cast<mpi::stencil_op<sten>&>(fop);
		grid_topo & topo = fopd.grid();
		MsgCtx *ctx = (MsgCtx*) fopd.halo_ctx;

		kc = topo.level() + 1;
		nog = topo.nlevel();
		kf = kc + 1;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		if (std::is_same<sten, seven_pt>::value) {
			MPI_BMG3_SymStd_SETUP_ITLI07_ex(kf, kc, fopd.data(), cop.data(), Pd.data(),
			                                fop.len(0), fop.len(1), fop.len(2),
			                                cop.len(0), cop.len(1), cop.len(2),
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
			MPI_BMG3_SymStd_SETUP_ITLI27_ex(kf, kc, fopd.data(), cop.data(), Pd.data(),
			                                fop.len(0), fop.len(1), fop.len(2),
			                                cop.len(0), cop.len(1), cop.len(2),
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

#endif
