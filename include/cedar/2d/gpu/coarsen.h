#ifndef CEDAR_2D_GPU_COARSEN_H
#define CEDAR_2D_GPU_COARSEN_H

#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>
#include <cedar/kernels/coarsen_op.h>
#include <cedar/2d/gpu/types.h>

using real_t = cedar::real_t;
using len_t = cedar::len_t;
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_ITLI_ex.f90.hpp>

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_ITLI_ex(int kf, int kc, real_t *SO, real_t *SOC, real_t *CI,
	                                   len_t IIF, len_t JJF, len_t IIC, len_t JJC, len_t iGs, len_t jGs,
	                                   int nog, int ifd, int nstencil,
	                                   void *halof);
}

namespace cedar { namespace cdr2 { namespace gpu { namespace mpi {

class galerkin : public kernels::coarsen_op<stypes>
{
	void run(const prolong_op & P,
	         const stencil_op<five_pt> & fop,
	         stencil_op<nine_pt> & cop) override
	{
		this->run_impl(P, fop, cop);
	}
	void run(const prolong_op & P,
	         const stencil_op<nine_pt> & fop,
	         stencil_op<nine_pt> & cop) override
	{
		this->run_impl(P, fop, cop);
	}


	template<class sten>
	void run_impl(const prolong_op & P,
	              const stencil_op<sten> & fop,
	              stencil_op<nine_pt> & cop)
	{
		int ifd, nstencil;
		int kf, kc, nog;
		auto & Pd = const_cast<prolong_op&>(P);
		auto & fopd = const_cast<mpi::stencil_op<sten>&>(fop);
		grid_topo & topo = fopd.grid();

		nstencil = stencil_ndirs<sten>::value;
		if (std::is_same<five_pt, sten>::value)
			ifd = 1;
		else
			ifd = 0;

		kc = topo.level();
		nog = topo.nlevel();
		kf = kc + 1;

                auto fopdb = fopd.to_buffer();

                Pd.ensure_gpu();
                fopd.ensure_gpu();
                cop.ensure_gpu();

                void* halof = services->fortran_handle<halo_exchange>();

                MPI_BMG2_SymStd_SETUP_ITLI_ex<ftl::device::GPU>(
                    kf, kc, fopd, cop, Pd, fop.len(0), fop.len(1), cop.len(0), cop.len(1),
                    topo.is(0), topo.is(1), nog, ifd, nstencil, halof);

                    // MPI_BMG2_SymStd_SETUP_ITLI_ex(kf, kc, fopd.data(), cop.data(), Pd.data(),
                    //                               fop.len(0), fop.len(1), cop.len(0), cop.len(1),
                    //                               topo.is(0), topo.is(1),
                    //                               nog, ifd, nstencil,
                    //                               services->fortran_handle<halo_exchange>());

                    cop.mark_cpu_dirty(true);

                // std::cerr << "coarse operator" << std::endl;
                // std::cerr << "has cpu: " << cop.has_cpu() << std::endl;
                // std::cerr << "has gpu: " << cop.has_gpu() << std::endl;
                // std::cerr << "cpu ptr: " << cop.to_flat_buffer().get_host_impl()->get_host_pointer() << std::endl;
                // std::cerr << "dev ptr: " << cop.to_flat_buffer().get_dev_impl().get() << std::endl;

                auto fopb = fopd.to_buffer();
                fopb.reshape({ fopb.len(0) + 1, fopb.len(1) + 1, fopb.len(2) });

                auto copb = cop.to_buffer();
                copb.reshape({ copb.len(0) + 1, copb.len(1) + 1, copb.len(2) });

                std::cerr << " == Galerkin product == " << std::endl;
                std::cerr << "Fine operator: " << std::endl << fopb << std::endl;
                std::cerr << "Fine ptr: " << fopb.data() << std::endl;
                std::cerr << "Fine obj ptr: " << (&fop) << std::endl;
                std::cerr << "Coarse operator: " << std::endl << copb << std::endl;
                std::cerr << "Coarse ptr: " << copb.data() << std::endl;
                std::cerr << "Coarse obj ptr: " << (&cop) << std::endl;
                std::cerr << " ====================== " << std::endl;
	}

};

}}}}

#endif
