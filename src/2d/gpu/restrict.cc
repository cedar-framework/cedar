#include <cedar/2d/ftn/BMG_parameters_c.h>
#include <cedar/2d/gpu/restrict.h>

#include <typeinfo>

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

namespace cedar::cdr2::gpu::mpi {
    void restrict_f90::update_periodic(mpi::grid_func & q,
                                       const grid_topo & topof,
                                       const grid_topo & topoc,
                                       std::array<bool, 3> periodic) {
        std::array<len_t, 2> ngf({topof.nglobal(0), topof.nglobal(1)});
        std::array<len_t, 2> ngc({topoc.nglobal(0), topoc.nglobal(1)});

        q.ensure_cpu();

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

        q.mark_cpu_dirty(true);
    }


    void restrict_f90::run(const restrict_op & R,
                           const grid_func & fine,
                           grid_func & coarse) {
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

        fined.ensure_gpu();
        coarse.ensure_gpu();
        P.ensure_gpu();

        // std::cerr << "coarse residual" << std::endl;
        // std::cerr << "has cpu: " << coarse.has_cpu() << std::endl;
        // std::cerr << "has gpu: " << coarse.has_gpu() << std::endl;
        // std::cerr << "cpu ptr: " << coarse.to_flat_buffer().get_host_impl()->get_host_pointer() << std::endl;
        // std::cerr << "dev ptr: " << coarse.to_flat_buffer().get_dev_impl().get() << std::endl;

        auto coarseb = coarse.to_buffer();
        std::cerr << " == Restriction == " << std::endl;
        std::cerr << "Coarse residual before restriction: " << std::endl << coarseb << std::endl;
        std::cerr << coarse.to_flat_buffer().tostr_unformatted() << std::endl;

        std::cerr << "Host implementation type: " << typeid(*coarse.to_flat_buffer().get_host_impl().get()).name() << std::endl;

        MPI_BMG2_SymStd_restrict<ftl::device::GPU>(
            kf, kc, nog, fined, coarse, P,
            fined.len(0), fined.len(1),
            coarse.len(0), coarse.len(1),
            fine_topo.is(0), fine_topo.is(1));

        // MPI_BMG2_SymStd_restrict(kf, kc, nog,
        //                          fined.data(), coarse.data(), P.data(),
        //                          fined.len(0), fined.len(1),
        //                          coarse.len(0), coarse.len(1),
        //                          fine_topo.is(0), fine_topo.is(1));

        // coarse.mark_cpu_dirty(true);

        auto fineb = fined.to_buffer();
        auto Pb = P.to_buffer();

        std::cerr << "Residual: " << std::endl << fineb << std::endl;
        std::cerr << "Restricted residual: " << std::endl << coarseb << std::endl;
        std::cerr << coarse.to_flat_buffer().tostr_unformatted() << std::endl;
        std::cerr << coarse.data() << std::endl;
        std::cerr << " ================= " << std::endl;
    }

}
