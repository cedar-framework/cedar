#include <boxmg/2d/ftn/BMG_parameters_c.h>
#include <boxmg/3d/mpi/halo.h>

#include <boxmg/3d/inter/restrict.h>

extern "C" {
	using namespace boxmg;
	void MPI_BMG3_SymStd_restrict(int kfg, int kcg,
	                              real_t *q, real_t *qc, real_t *ci,
	                              len_t nx, len_t ny, len_t nz,
	                              len_t nxc, len_t nyc, len_t nzc,
	                              len_t igs, len_t jgs, len_t kgs);
}

namespace boxmg { namespace bmg3 { namespace kernel {

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
		MPI_BMG3_SymStd_restrict(kf, kc,
		                         fined.data(), coarse.data(), P.data(),
		                         fined.len(0), fined.len(1), fined.len(2),
		                         coarse.len(0), coarse.len(1), coarse.len(2),
		                         topo.is(0), topo.is(1), topo.is(2));
	}
}

}}}
