#include "cedar/2d/ftn/BMG_parameters_c.h"
#include "cedar/2d/inter/mpi/prolong_op.h"

#include "cedar/2d/mpi/halo.h"
#include "cedar/2d/inter/restrict.h"

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_restrict(int kf, int kc, int nog,
	                              real_t *Q, real_t *QC, real_t *CI,
	                              len_t nx, len_t ny, len_t nxc, len_t nyc,
	                              len_t iGs, len_t jGs,
	                              len_t *iWork, len_t NMSGi, int *pMSG,
	                              real_t *msg_buffer, len_t NMSGr, int MPICOMM);
}

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	void mpi_fortran_restrict(const inter::mpi::restrict_op & R,
	                          const mpi::grid_func & fine,
	                          mpi::grid_func & coarse)
	{
		int kf, kc, nog;
		auto & Rd = const_cast<inter::mpi::restrict_op&>(R);
		inter::mpi::prolong_op & P = Rd.getP();
		grid_topo & topo = P.grid();
		const grid_topo & fine_topo = fine.grid();
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
