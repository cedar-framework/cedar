#include "cedar/2d/ftn/BMG_parameters_c.h"
#include "cedar/2d/inter/mpi/prolong_op.h"

#include "cedar/2d/inter/restrict.h"

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_restrict(int kf, int kc, int nog,
	                              real_t *Q, real_t *QC, real_t *CI,
	                              len_t nx, len_t ny, len_t nxc, len_t nyc,
	                              len_t iGs, len_t jGs);
}


namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	static void update_periodic(mpi::grid_func & q,
	                            const grid_topo & topof,
	                            const grid_topo & topoc,
	                            std::array<bool, 3> periodic)
	{
		std::array<len_t, 2> ngf({topof.nglobal(0), topof.nglobal(1)});
		std::array<len_t, 2> ngc({topoc.nglobal(0), topoc.nglobal(1)});

		if (periodic[0] and ((ngf[0]/2 + 1) != ngc[0])) {
			if (topof.coord(0) == 0) {
				for (auto j : q.grange(1)) {
					q(0,j) = 0;
				}
			}
			if (topof.coord(0) == (topof.nproc(0) - 1)) {
				for (auto j : q.grange(1)) {
					q(q.len(0)-1,j) = 0;
				}
			}
		}

		if (periodic[1] and ((ngf[1]/2 + 1) != ngc[1])) {
			if (topof.coord(1) == 0) {
				for (auto i : q.grange(0)) {
					q(i,0) = 0;
				}
			}
			if (topof.coord(1) == (topof.nproc(1) - 1)) {
				for (auto i : q.grange(0)) {
					q(i,q.len(1)-1) = 0;
				}
			}
		}
	}


	void mpi_fortran_restrict(const kernel_params & params,
	                          const inter::mpi::restrict_op & R,
	                          const mpi::grid_func & fine,
	                          mpi::grid_func & coarse)
	{
		int kf, kc, nog;
		auto & Rd = const_cast<inter::mpi::restrict_op&>(R);
		inter::mpi::prolong_op & P = Rd.getP();
		grid_topo & topo = P.grid();
		const grid_topo & fine_topo = fine.grid();
		const grid_topo & coarse_topo = coarse.grid();
		auto & fined = const_cast<mpi::grid_func&>(fine);

		nog = kf = topo.nlevel();
		kc = topo.level();

		// conditionally zero periodic entries
		update_periodic(fined, fine_topo, coarse_topo, params.periodic);

		MPI_BMG2_SymStd_restrict(kf, kc, nog,
		                         fined.data(), coarse.data(), P.data(),
		                         fined.len(0), fined.len(1),
		                         coarse.len(0), coarse.len(1),
		                         fine_topo.is(0), fine_topo.is(1));
	}
}

}}}
