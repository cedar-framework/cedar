#ifndef CEDAR_2D_GPU_INTERP_H
#define CEDAR_2D_GPU_INTERP_H

#include <cedar/2d/gpu/types.h>
#include <cedar/kernels/interp_add.h>
#include <cedar/kernels/setup_interp.h>
#include <cedar/2d/gpu/kernel_manager.h>

using real_t = cedar::real_t;
using len_t = cedar::len_t;
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_interp_OI.f90.hpp>

extern "C" {
    void MPI_BMG2_SymStd_SETUP_interp_OI(int kf, int kc, real_t *so, real_t *ci,
                                         len_t iif, len_t jjf, len_t iic, len_t jjc,
                                         int nog, int nogm, len_t *igrd,
                                         int ifd, int nstncl, int ibc,
                                         void *halof);

    void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr2 { namespace gpu { namespace mpi {

class interp_f90 : public kernels::interp_add<stypes>
{
	void run(const prolong_op & P,
	         const grid_func & coarse,
	         const grid_func & residual,
	         grid_func & fine) override;

};


class setup_interp_f90 : public kernels::setup_interp<stypes>
{
	void run(const stencil_op<five_pt> & fop,
	         const stencil_op<nine_pt> & cop,
	         prolong_op & P) override { this->run_impl(fop, cop, P); }
	void run(const stencil_op<nine_pt> & fop,
	         const stencil_op<nine_pt> & cop,
	         prolong_op & P) override { this->run_impl(fop, cop, P); }

	template<class sten>
	void run_impl(const stencil_op<sten> & fop,
	              const stencil_op<nine_pt> & cop,
	              prolong_op & P)
	{
		int ifd, nstencil, jpn;
		int kf, kc, nog;

		auto & fopd = const_cast<stencil_op<sten>&>(fop);
		grid_topo & topo = fopd.grid();

		store_fine_op(fopd, P);

		if (stencil_ndirs<sten>::value == 3) {
			ifd = 1;
			nstencil = 3;
		} else {
			ifd = 0;
			nstencil = 5;
		}

		kc = topo.level();
		nog = topo.nlevel();
		kf = kc + 1;

		BMG_get_bc(params->per_mask(), &jpn);

                void* halof = services->fortran_handle<halo_exchange>();

                auto igrd_vec = topo.get_igrd();
                ftl::FlatBuffer<len_t, len_t> igrd(igrd_vec->data(), igrd_vec->size());

                fopd.ensure_gpu();
                const_cast<stencil_op<nine_pt>&>(cop).ensure_gpu();
                P.ensure_gpu();

		MPI_BMG2_SymStd_SETUP_interp_OI<ftl::device::GPU>(
                    kf, kc, fopd, P,
                    fop.len(0), fop.len(1), cop.len(0), cop.len(1),
                    nog, nog, igrd, ifd, nstencil, jpn, halof);

                // std::cerr << "interpolation operator" << std::endl;
                // std::cerr << "has cpu: " << P.has_cpu() << std::endl;
                // std::cerr << "has gpu: " << P.has_gpu() << std::endl;
                // std::cerr << "cpu ptr: " << P.to_flat_buffer().get_host_impl()->get_host_pointer() << std::endl;
                // std::cerr << "dev ptr: " << P.to_flat_buffer().get_dev_impl().get() << std::endl;

                auto Pb = P.to_buffer();
                std::cerr << " == Form Interpolation == " << std::endl;
                std::cerr << "Operator: " << std::endl << Pb << std::endl;
                std::cerr << " =================== " << std::endl;
	}


	inline void store_fine_op(stencil_op<five_pt> & fop,
	                          prolong_op & P)
	{
		P.fine_op_five = &fop;
		P.fine_is_five = true;
	}

	inline void store_fine_op(stencil_op<nine_pt> & fop,
		                          prolong_op & P)
	{
		P.fine_op_nine = &fop;
		P.fine_is_five = false;
	}
};

}}}}

#endif
