#include <cedar/2d/ftn/BMG_parameters_c.h>
#include <cedar/2d/gpu/interp.h>

using real_t = cedar::real_t;
using len_t = cedar::len_t;

#include <src/2d/ftn/mpi/BMG2_SymStd_interp_add.f90.hpp>

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_interp_add(len_t kc, len_t kf, len_t nog,
	                                real_t *Q, real_t *QC, real_t *RES, real_t *SO,
	                                int nstncl, real_t *CI, len_t iic, len_t jjc,
	                                len_t iif, len_t jjf, len_t iGs, len_t jGs,
	                                void *halof);
}

namespace cedar::cdr2::gpu::mpi {

    void interp_f90::run(const prolong_op & P,
                         const grid_func & coarse,
                         const mpi::grid_func & residual,
                         grid_func & fine) {
        int nstencil, kf, kc, nog;

        auto & Pd = const_cast<prolong_op&>(P);
        grid_func & coarsed = const_cast<grid_func&>(coarse);
        grid_func & res = const_cast<grid_func&>(residual);
        grid_topo & topo = Pd.grid();

        topo_ptr topof;

        ftl::Buffer<real_t, len_t> fop_data;

        if (Pd.fine_is_five) {
            nstencil = 3;
            fop_data = Pd.fine_op_five->to_buffer();
            topof = Pd.fine_op_five->grid_ptr();
        } else {
            nstencil = 5;
            fop_data = Pd.fine_op_nine->to_buffer();
            topof = Pd.fine_op_nine->grid_ptr();
        }

        kc = topo.level() + 1;
        nog = topo.nlevel();
        kf = kc + 1;

        void* halof = services->fortran_handle<halo_exchange>();


        // Pd.ensure_cpu();
        // coarsed.ensure_cpu();
        // res.ensure_cpu();
        // fine.ensure_cpu();

        std::cerr << "Residual dirty: (host " << res.is_cpu_dirty() <<
            "   dev " << res.is_gpu_dirty() << ")" << std::endl;

        Pd.ensure_cpu();
        coarsed.ensure_cpu();
        //res.ensure_cpu();
        fine.ensure_cpu();

        // MPI_BMG2_SymStd_interp_add(kc, kf, nog,
        //                            fine.data(), coarsed.data(), res.data(),
        //                            fop_data.data(), nstencil,
        //                            Pd.data(),
        //                            coarsed.len(0), coarsed.len(1),
        //                            fine.len(0), fine.len(1),
        //                            topof->is(0), topof->is(1),
        //                            services->template fortran_handle<halo_exchange>());

        MPI_BMG2_SymStd_interp_add<ftl::device::GPU>(
            kc, kf, nog, fine, coarsed, res,
            fop_data, nstencil, Pd,
            coarsed.len(0), coarsed.len(1),
            fine.len(0), fine.len(1),
            topof->is(0), topof->is(1), halof);

        fine.mark_cpu_dirty(true);

        auto coarseb = coarsed.to_buffer();
        auto fineb = fine.to_buffer();
        auto resb = res.to_buffer();

        std::cerr << " == Interpolation == " << std::endl;
        std::cerr << "Coarse solution: " << std::endl << coarseb << std::endl;
        std::cerr << "Residual : " << std::endl << resb << std::endl;
        std::cerr << "Fine solution: " << std::endl << fineb << std::endl;
        std::cerr << " =================== " << std::endl;
    }

}
