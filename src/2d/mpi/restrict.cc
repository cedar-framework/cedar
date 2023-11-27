#include <cedar/2d/ftn/BMG_parameters_c.h>
#include <cedar/2d/mpi/restrict.h>

using real_t = cedar::real_t;
using len_t = cedar::len_t;

#include <src/2d/ftn/mpi/BMG2_SymStd_restrict.f90.hpp>

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_restrict(int kf, int kc, int nog,
	                              real_t *Q, real_t *QC, real_t *CI,
	                              len_t nx, len_t ny, len_t nxc, len_t nyc,
	                              len_t iGs, len_t jGs);
}

namespace cedar { namespace cdr2 { namespace mpi {

        template <typename device>
	void restrict_f90<device>::update_periodic(mpi::grid_func & q,
	                            const grid_topo & topof,
	                            const grid_topo & topoc,
	                            std::array<bool, 3> periodic)
	{
		std::array<len_t, 2> ngf({topof.nglobal(0), topof.nglobal(1)});
		std::array<len_t, 2> ngc({topoc.nglobal(0), topoc.nglobal(1)});

		if (periodic[0] and ((ngf[0]/2 + 1) != ngc[0])) {
                    std::cerr << "Ensuring CPU for periodic restriction" << std::endl;
                    q.ensure_cpu();
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
                        q.mark_cpu_dirty(true);
		}

		if (periodic[1] and ((ngf[1]/2 + 1) != ngc[1])) {
                    std::cerr << "Ensuring CPU for periodic restriction" << std::endl;
                    q.ensure_cpu();
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
                        q.mark_cpu_dirty(true);
		}
	}

        template
        void restrict_f90<cedar::cpu>::update_periodic(
            mpi::grid_func &, const grid_topo &, const grid_topo &, std::array<bool, 3>);

        template
        void restrict_f90<cedar::gpu>::update_periodic(
            mpi::grid_func &, const grid_topo &, const grid_topo &, std::array<bool, 3>);

        template <typename device>
	void restrict_f90<device>::run(const restrict_op & R,
                       const grid_func & fine,
                       grid_func & coarse)
	{
		int kf, kc, nog;
		auto & Rd = const_cast<restrict_op&>(R);
		prolong_op & P = Rd.getP();
		grid_topo & topo = P.grid();
		const grid_topo & fine_topo = fine.grid();
		const grid_topo & coarse_topo = coarse.grid();
		auto & fined = const_cast<mpi::grid_func&>(fine);

		nog = kf = topo.nlevel();
		kc = topo.level();

		// conditionally zero periodic entries
		update_periodic(fined, fine_topo, coarse_topo, params->periodic);

                fined.template ensure<device>();
                coarse.template ensure<device>();
                P.template ensure<device>();

                if (device::is_gpu()) {
                    MPI_BMG2_SymStd_restrict<cedar::gpu>(
                        kf, kc, nog, fined, coarse, P,
                        fined.len(0), fined.len(1),
                        coarse.len(0), coarse.len(1),
                        fine_topo.is(0), fine_topo.is(1));
                } else {
                    MPI_BMG2_SymStd_restrict(kf, kc, nog,
                                             fined.data(), coarse.data(), P.data(),
                                             fined.len(0), fined.len(1),
                                             coarse.len(0), coarse.len(1),
                                             fine_topo.is(0), fine_topo.is(1));
                    coarse.template mark_dirty<device>();
                }
	}

        template
        void restrict_f90<cedar::cpu>::run(
            const restrict_op &, const grid_func &, grid_func &);

        template
        void restrict_f90<cedar::gpu>::run(
            const restrict_op &, const grid_func &, grid_func &);
}}}
