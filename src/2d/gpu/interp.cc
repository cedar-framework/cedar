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

        ftl::Buffer<real_t, len_t> fop_data;
        topo_ptr topof;

        if (Pd.fine_is_five) {
            nstencil = 3;
            fop_data = Pd.fine_op_five->to_buffer();
            topof = Pd.fine_op_five->grid_ptr();
            std::cerr << "Fine op object ptr: " << Pd.fine_op_five << std::endl;
            std::cerr << "Fine op data ptr: " << Pd.fine_op_five->data() << std::endl;
        } else {
            nstencil = 5;
            fop_data = Pd.fine_op_nine->to_buffer();
            topof = Pd.fine_op_nine->grid_ptr();
            std::cerr << "Fine op object ptr: " << Pd.fine_op_nine << std::endl;
            std::cerr << "Fine op data ptr: " << Pd.fine_op_nine->data() << std::endl;
        }

        kc = topo.level() + 1;
        nog = topo.nlevel();
        kf = kc + 1;

        void* halof = services->fortran_handle<halo_exchange>();

        Pd.ensure_cpu();
        coarsed.ensure_cpu();
        fine.ensure_cpu();

        MPI_BMG2_SymStd_interp_add<ftl::device::GPU>(
            kc, kf, nog, fine, coarsed, res,
            fop_data, nstencil, Pd,
            coarsed.len(0), coarsed.len(1),
            fine.len(0), fine.len(1),
            topof->is(0), topof->is(1), halof);

        // MPI_BMG2_SymStd_interp_add(kc, kf, nog,
        //                            fine.data(), coarsed.data(), res.data(),
        //                            fop_data.data(), nstencil,
        //                            Pd.data(),
        //                            coarsed.len(0), coarsed.len(1),
        //                            fine.len(0), fine.len(1),
        //                            topof->is(0), topof->is(1),
        //                            services->template fortran_handle<halo_exchange>());

        fine.mark_cpu_dirty(true);

        auto coarseb = coarsed.to_buffer();
        auto fineb = fine.to_buffer();
        auto resb = res.to_buffer();
        auto Pb = Pd.to_buffer();
        ftl::Buffer<real_t, len_t> fopb = fop_data;
        fopb.reshape({ fopb.len(0) + 1, fopb.len(1) + 1, fopb.len(2) });

        // std::cerr << "fine-grid solution" << std::endl;
        // std::cerr << "has cpu: " << fine.has_cpu() << std::endl;
        // std::cerr << "has gpu: " << fine.has_gpu() << std::endl;
        // std::cerr << "cpu ptr: " << fine.to_flat_buffer().get_host_impl()->get_host_pointer() << std::endl;
        // std::cerr << "dev ptr: " << fine.to_flat_buffer().get_dev_impl().get() << std::endl;

        std::cerr << " == Interpolation == " << std::endl;
        std::cerr << "Fine operator: " << std::endl << fopb << std::endl;
        std::cerr << "Fine ptr: " << fopb.data() << std::endl;
        std::cerr << "Coarse solution: " << std::endl << coarseb << std::endl;
        std::cerr << "Residual : " << std::endl << resb << std::endl;
        std::cerr << "Fine solution: " << std::endl << fineb << std::endl;
        std::cerr << " =================== " << std::endl;
    }

}
