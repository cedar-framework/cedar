#include "boxmg/2d/ftn/BMG_parameters_c.h"
#include "boxmg/2d/inter/mpi/prolong_op.h"

#include "boxmg/2d/mpi/halo.h"
#include "boxmg/2d/inter/restrict.h"

extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_restrict(real_t*, real_t*, real_t*,
	                          int, int, int, int, int);
	void MPI_BMG2_SymStd_restrict(int kf, int kc, int nog,
	                              real_t *Q, real_t *QC, real_t *CI,
	                              len_t nx, len_t ny, len_t nxc, len_t nyc,
	                              len_t iGs, len_t jGs,
	                              len_t *iWork, len_t NMSGi, int *pMSG,
	                              real_t *msg_buffer, len_t NMSGr, int MPICOMM);
}

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void fortran_restrict(const inter::restrict_op & R,
	                      const grid_func & fine,
	                      grid_func & coarse)
	{
		using namespace boxmg::bmg2d;
		int ibc;

		auto & fined = const_cast<grid_func&>(fine);
		auto & Rd = const_cast<inter::restrict_op&>(R);
		inter::prolong_op & P = Rd.getP();
		ibc = BMG_BCs_definite;

		BMG2_SymStd_restrict(fined.data(), coarse.data(),
		                     P.data(), fined.len(0), fined.len(1),
		                     coarse.len(0), coarse.len(1), ibc);
	}

	void mpi_fortran_restrict(const inter::mpi::restrict_op & R,
	                          const mpi::grid_func & fine,
	                          mpi::grid_func & coarse)
	{
		int kf, kc, nog;
		auto & Rd = const_cast<inter::mpi::restrict_op&>(R);
		inter::mpi::prolong_op & P = Rd.getP();
		mpi::grid_topo & topo = P.grid();
		const mpi::grid_topo & fine_topo = fine.grid();
		MsgCtx *ctx = (MsgCtx*) P.halo_ctx;
		auto & fined = const_cast<mpi::grid_func&>(fine);

		nog = kf = topo.nlevel();
		kc = topo.level();

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		MPI_BMG2_SymStd_restrict(kf, kc, nog,
		                         fined.data(), coarse.data(), P.data(),
		                         fined.len(0), fined.len(1),
		                         coarse.len(0), coarse.len(1),
		                         fine_topo.is(0), fine_topo.is(1),
		                         ctx->msg_geom.data(), ctx->msg_geom.size(),
		                         ctx->pMSG.data(), ctx->msg_buffer.data(),
		                         ctx->msg_buffer.size(), fcomm);
	}
}

}}}
