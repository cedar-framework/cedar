#include <cedar/2d/ftn/BMG_parameters_c.h>

#include <cedar/3d/mpi/restrict.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_restrict(int kfg, int kcg,
	                              real_t *q, real_t *qc, real_t *ci,
	                              len_t nx, len_t ny, len_t nz,
	                              len_t nxc, len_t nyc, len_t nzc,
	                              len_t igs, len_t jgs, len_t kgs);
}

namespace cedar { namespace cdr3 { namespace mpi {

void restrict_f90::run(const restrict_op & R,
                       const grid_func & fine,
                       grid_func & coarse)
{
	int kf, kc;
	auto & Rd = const_cast<restrict_op&>(R);
	prolong_op & P = Rd.getP();
	grid_topo & topo = P.grid();
	const grid_topo & fine_topo = fine.grid();
	auto & fined = const_cast<grid_func&>(fine);

	kf = topo.nlevel();
	kc = topo.level();

	MPI_BMG3_SymStd_restrict(kf, kc,
	                         fined.data(), coarse.data(), P.data(),
	                         fined.len(0), fined.len(1), fined.len(2),
	                         coarse.len(0), coarse.len(1), coarse.len(2),
	                         fine_topo.is(0), fine_topo.is(1), fine_topo.is(2));
}


}}}
