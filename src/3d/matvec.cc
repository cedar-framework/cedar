#include <cedar/3d/matvec.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/halo.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_UTILS_matvec(int kg, real_t *so, real_t *qf,
	                                  real_t *q, len_t ii, len_t jj, len_t kk,
	                                  int nog, int ifd, int nstncl,
	                                  len_t *iwork, int *pMSG, real_t *msg_buffer,
	                                  int mpicomm);
}

namespace cedar { namespace cdr3 { namespace kernel {
namespace impls
{
	namespace mpi = cedar::cdr3::mpi;

	void matvec(const mpi::stencil_op & so,
	            const mpi::grid_func & x,
	            mpi::grid_func & b)
	{
		using namespace cedar::cdr3;
		int kg, ifd, nstencil, nog;

		mpi::stencil_op & sod = const_cast<mpi::stencil_op&>(so);
		mpi::grid_func & xd = const_cast<mpi::grid_func&>(x);
		grid_topo & topo = sod.grid();
		MsgCtx *ctx = (MsgCtx*) sod.halo_ctx;
		grid_stencil & sten = sod.stencil();

		nog = topo.nlevel();
		kg = topo.level()+1;
		if (sten.five_pt()) {
			ifd = 1;
			nstencil = 3;
		} else {
			ifd = 0;
			nstencil = 5;
		}

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		MPI_BMG3_SymStd_UTILS_matvec(kg, sod.data(), b.data(), xd.data(),
		                             sten.len(0), sten.len(1), sten.len(2),
		                             nog, ifd, nstencil,
		                             ctx->msg_geom.data(), ctx->pMSG.data(),
		                             ctx->msg_buffer.data(), fcomm);
	}
}

}}}
