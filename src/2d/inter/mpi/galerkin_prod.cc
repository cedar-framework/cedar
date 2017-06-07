#include "cedar/2d/ftn/BMG_parameters_c.h"
#include "cedar/2d/mpi/stencil_op.h"
#include "cedar/2d/inter/mpi/prolong_op.h"

#include "cedar/2d/mpi/halo.h"
#include "cedar/2d/inter/galerkin_prod.h"


extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_ITLI_ex(int kf, int kc, real_t *SO, real_t *SOC, real_t *CI,
	                                   len_t IIF, len_t JJF, len_t IIC, len_t JJC, len_t iGs, len_t jGs,
	                                   int nog, int ifd, int nstencil,
	                                   len_t *iWork, len_t NMSGi, int *pMSGSO,
	                                   real_t *msg_buffer, len_t NMSGr, int MPICOMM);
}


namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr2;

	template<>
	void mpi_galerkin_prod(const kernel_params & params,
	                       const inter::mpi::prolong_op & P,
	                       const mpi::stencil_op<five_pt> & fop,
	                       mpi::stencil_op<nine_pt> & cop)
	{
		int ifd, nstencil;
		int kf, kc, nog;
		auto & Pd = const_cast<inter::mpi::prolong_op&>(P);
		auto & fopd = const_cast<mpi::stencil_op<five_pt>&>(fop);
		grid_topo & topo = fopd.grid();
		MsgCtx *ctx = (MsgCtx*) fopd.halo_ctx;

		ifd = 1;
		nstencil = 3;

		kc = topo.level();
		nog = topo.nlevel();
		kf = kc + 1;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_BMG2_SymStd_SETUP_ITLI_ex(kf, kc, fopd.data(), cop.data(), Pd.data(),
		                              fop.len(0), fop.len(1), cop.len(0), cop.len(1),
		                              topo.is(0), topo.is(1),
		                              nog, ifd, nstencil,
		                              ctx->msg_geom.data(), ctx->msg_geom.size(), ctx->pMSGSO.data(),
		                              ctx->msg_buffer.data(), ctx->msg_buffer.size(), fcomm);
	}


	template<>
	void mpi_galerkin_prod(const kernel_params & params,
	                       const inter::mpi::prolong_op & P,
	                       const mpi::stencil_op<nine_pt> & fop,
	                       mpi::stencil_op<nine_pt> & cop)
	{
		int ifd, nstencil;
		int kf, kc, nog;
		auto & Pd = const_cast<inter::mpi::prolong_op&>(P);
		auto & fopd = const_cast<mpi::stencil_op<nine_pt>&>(fop);
		grid_topo & topo = fopd.grid();
		MsgCtx *ctx = (MsgCtx*) fopd.halo_ctx;

		ifd = 0;
		nstencil = 5;

		kc = topo.level();
		nog = topo.nlevel();
		kf = kc + 1;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_BMG2_SymStd_SETUP_ITLI_ex(kf, kc, fopd.data(), cop.data(), Pd.data(),
		                              fop.len(0), fop.len(1), cop.len(0), cop.len(1),
		                              topo.is(0), topo.is(1),
		                              nog, ifd, nstencil,
		                              ctx->msg_geom.data(), ctx->msg_geom.size(), ctx->pMSGSO.data(),
		                              ctx->msg_buffer.data(), ctx->msg_buffer.size(), fcomm);
	}

}

}}}
