#include "fortran/BMG_parameters_c.h"
#include "core/mpi/stencil_op.h"
#include "inter/mpi/prolong_op.h"

#include "halo.h"
#include "galerkin_prod.h"


extern "C" {
	void bmg2_symstd_setup_itli_ex(int*,int*,double*,double*,double*,
	                               int*,int*,int*,int*, int*, int*, int*,int*);
	using namespace boxmg;
	void MPI_BMG2_SymStd_SETUP_ITLI_ex(int kf, int kc, real_t *SO, real_t *SOC, real_t *CI,
	                                   len_t IIF, len_t JJF, len_t IIC, len_t JJC, len_t iGs, len_t jGs,
	                                   int nog, int ifd, int nstencil,
	                                   len_t *iWork, len_t NMSGi, int *pMSGSO,
	                                   real_t *msg_buffer, len_t NMSGr, int MPICOMM);
}


namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	using namespace boxmg::bmg2d;

	void galerkin_prod(int kf, int kc, int nog,
	                   const inter::prolong_op & P,
	                   const stencil_op & fop,
	                   stencil_op & cop)
	{
		using namespace boxmg::bmg2d::core;
		int iif, jjf, iic, jjc, ifd;
		int nstencil, ipn;

		const grid_stencil &fsten = fop.stencil();
		grid_stencil & csten = cop.stencil();
		stencil_op &fopd = const_cast<stencil_op&>(fop);
		inter::prolong_op &Pd = const_cast<inter::prolong_op&>(P);

		iif = fsten.len(0);
		jjf = fsten.len(1);
		iic = csten.len(0);
		jjc = csten.len(1);

		if (fsten.five_pt()) {
			ifd = 0;
			nstencil = 5;
		} else {
			ifd = 1;
			nstencil = 9;
		}

		ipn = BMG_BCs_definite;

		bmg2_symstd_setup_itli_ex(&kf, &kc, fopd.data(), cop.data(), Pd.data(),
		                          &iif, &jjf, &iic, &jjc, &nog, &ifd, &nstencil, &ipn);

	}


	void mpi_galerkin_prod(int kf, int kc, int nog,
	                       const inter::prolong_op & P,
	                       const stencil_op & fop,
	                       stencil_op & cop)
	{
		int ifd, nstencil;
		inter::mpi::prolong_op & Pd = const_cast<inter::mpi::prolong_op&>(dynamic_cast<const inter::mpi::prolong_op&>(P));
		mpi::stencil_op & fopd = const_cast<mpi::stencil_op&>(dynamic_cast<const mpi::stencil_op&>(fop));
		mpi::stencil_op & copd = dynamic_cast<mpi::stencil_op&>(cop);
		mpi::grid_topo & topo = fopd.grid();
		grid_stencil & fsten = fopd.stencil();
		grid_stencil & csten = copd.stencil();
		MsgCtx *ctx = (MsgCtx*) fopd.halo_ctx;

		if (fsten.five_pt()) {
			ifd = 0;
			nstencil = 3;
		} else {
			ifd = 1;
			nstencil = 5;
		}

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_BMG2_SymStd_SETUP_ITLI_ex(kf, kc, fopd.data(), copd.data(), Pd.data(),
		                              fsten.len(0), fsten.len(1), csten.len(0), csten.len(1),
		                              topo.is(0), topo.is(1),
		                              nog, ifd, nstencil,
		                              ctx->msg_geom.data(), ctx->msg_geom.size(), ctx->pMSGSO.data(),
		                              ctx->msg_buffer.data(), ctx->msg_buffer.size(), fcomm);
	}


}

}}}
