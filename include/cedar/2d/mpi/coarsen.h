#ifndef CEDAR_2D_MPI_COARSEN_H
#define CEDAR_2D_MPI_COARSEN_H

#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>
#include <cedar/kernels/coarsen_op.h>
#include <cedar/2d/mpi/types.h>

using real_t = cedar::real_t;
using len_t = cedar::len_t;
#include <cedar/device.h>
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_ITLI_ex.f90.hpp>

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_ITLI_ex(int kf, int kc, real_t *SO, real_t *SOC, real_t *CI,
	                                   len_t IIF, len_t JJF, len_t IIC, len_t JJC, len_t iGs, len_t jGs,
	                                   int nog, int ifd, int nstencil,
	                                   void *halof);
}

namespace cedar { namespace cdr2 { namespace mpi {

template <typename device=cedar::cpu>
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

                void* halof = services->fortran_handle<halo_exchange>();

                Pd.template ensure<device>();
                fopd.template ensure<device>();
                cop.template ensure<device>();

                if (device::is_gpu()) {
                    MPI_BMG2_SymStd_SETUP_ITLI_ex<cedar::gpu>(
                        kf, kc, fopd, cop, Pd, fop.len(0), fop.len(1), cop.len(0), cop.len(1),
                        topo.is(0), topo.is(1), nog, ifd, nstencil, halof);
                } else {
                    MPI_BMG2_SymStd_SETUP_ITLI_ex(
                        kf, kc, fopd.data(), cop.data(), Pd.data(),
                        fop.len(0), fop.len(1), cop.len(0), cop.len(1),
                        topo.is(0), topo.is(1), nog, ifd, nstencil, halof);
                    cop.template mark_dirty<device>();
                }



                // cop.ensure_cpu();
                // auto copb = cop.to_buffer();
                // copb.dev_to_host();
                // auto copb = cop.to_buffer();
                // copb.reshape({
                //         copb.get_shape()[0] + 1,
                //         copb.get_shape()[1] + 1,
                //         copb.get_shape()[2]
                //     });
                // std::cerr << "Coarse operator: " << std::endl << copb << std::endl;
	}

};

}}}

#endif
