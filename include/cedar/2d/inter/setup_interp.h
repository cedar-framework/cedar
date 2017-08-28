#ifndef CEDAR_2D_KERNEL_SETUP_INTERP_H
#define CEDAR_2D_KERNEL_SETUP_INTERP_H

#include <cedar/kernel_params.h>
#include <cedar/halo_exchanger.h>
#include "cedar/2d/grid_func.h"
#include "cedar/2d/stencil_op.h"
#include "cedar/2d/inter/prolong_op.h"
#include "cedar/2d/mpi/stencil_op.h"
#include "cedar/2d/inter/mpi/prolong_op.h"
#include "cedar/2d/mpi/halo.h"


extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_interp_OI(int kf, int kc, real_t *so, real_t *ci,
	                                     len_t iif, len_t jjf, len_t iic, len_t jjc,
	                                     int nog, int nogm, len_t *igrd,
	                                     int ifd, int nstncl, int ibc,
	                                     void *halof);
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	template <class sten>
	void setup_interp(const kernel_params & params,
	                  const stencil_op<sten> & fop,
	                  const stencil_op<nine_pt> & cop,
	                  inter::prolong_op & P);

	namespace mpi = cedar::cdr2::mpi;

	template<class sten>
		void store_fine_op(mpi::stencil_op<sten> & fop,
		                          inter::mpi::prolong_op & P);

	template<>
		inline void store_fine_op(mpi::stencil_op<five_pt> & fop,
		                                 inter::mpi::prolong_op & P)
	{
		P.fine_op_five = &fop;
		P.fine_is_five = true;
	}


	template<>
		inline void store_fine_op(mpi::stencil_op<nine_pt> & fop,
		                                 inter::mpi::prolong_op & P)
	{
		P.fine_op_nine = &fop;
		P.fine_is_five = false;
	}


	template <class sten>
	void mpi_setup_interp(const kernel_params & params,
	                      halo_exchanger *halof,
	                      const mpi::stencil_op<sten> & fop,
	                      const mpi::stencil_op<nine_pt> & cop,
	                      inter::mpi::prolong_op & P)
	{
		int ifd, nstencil, jpn;
		int kf, kc, nog;

		auto & fopd = const_cast<mpi::stencil_op<sten>&>(fop);
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

		BMG_get_bc(params.per_mask(), &jpn);

		MPI_BMG2_SymStd_SETUP_interp_OI(kf, kc, fopd.data(), P.data(),
		                                fop.len(0), fop.len(1), cop.len(0), cop.len(1),
		                                nog, nog, topo.IGRD(), ifd, nstencil, jpn, halof);
	}
}

}}}

#endif
