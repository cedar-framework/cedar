#ifndef CEDAR_2D_KERNEL_MPI_GALERKIN_PROD_H
#define CEDAR_2D_KERNEL_MPI_GALERKIN_PROD_H

#include <cedar/kernel_params.h>
#include <cedar/halo_exchanger.h>
#include "cedar/2d/stencil_op.h"
#include "cedar/2d/inter/prolong_op.h"
#include "cedar/2d/grid_func.h"
#include "cedar/2d/mpi/stencil_op.h"
#include "cedar/2d/inter/mpi/prolong_op.h"
#include "cedar/2d/ftn/BMG_parameters_c.h"


extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_ITLI_ex(int kf, int kc, real_t *SO, real_t *SOC, real_t *CI,
	                                   len_t IIF, len_t JJF, len_t IIC, len_t JJC, len_t iGs, len_t jGs,
	                                   int nog, int ifd, int nstencil,
	                                   void *halof);
}

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr2::mpi;
	template <class sten>
	void mpi_galerkin_prod(const kernel_params & params,
	                       halo_exchanger_base *halof,
	                       const inter::mpi::prolong_op & P,
	                       const mpi::stencil_op<sten> & fop,
	                       mpi::stencil_op<nine_pt> & cop)
	{
		int ifd, nstencil;
		int kf, kc, nog;
		auto & Pd = const_cast<inter::mpi::prolong_op&>(P);
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

		MPI_BMG2_SymStd_SETUP_ITLI_ex(kf, kc, fopd.data(), cop.data(), Pd.data(),
		                              fop.len(0), fop.len(1), cop.len(0), cop.len(1),
		                              topo.is(0), topo.is(1),
		                              nog, ifd, nstencil,
		                              halof);
	}
}

}}}

#endif
