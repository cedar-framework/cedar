#ifndef CEDAR_2D_KERNEL_MATVEC_H
#define CEDAR_2D_KERNEL_MATVEC_H

#include <cedar/kernel_params.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/halo_exchanger_base.h>


extern "C" {
	using namespace cedar;
	void BMG2_SymStd_UTILS_matvec(int k, real_t *SO, real_t *QF,
	                              real_t *Q, len_t II, len_t JJ,
	                              int kf, int ifd, int nstencil);
}

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr2::mpi;
	template <class sten>
	void matvec(const kernel_params & params,
	            halo_exchanger_base *halof,
	            const mpi::stencil_op<sten> & so,
	            const mpi::grid_func & x,
	            mpi::grid_func & y)
	{
		using namespace cedar::cdr2;
		int k, kf, ifd;
		int nstencil;

		auto & sod = const_cast<mpi::stencil_op<sten>&>(so);
		mpi::grid_func & xd = const_cast<mpi::grid_func&>(x);
		grid_topo & topo = sod.grid();

		k = topo.level()+1;
		kf = topo.nlevel();

		nstencil = stencil_ndirs<sten>::value;
		if (std::is_same<five_pt, sten>::value)
			ifd = 1;
		else
			ifd = 0;

		BMG2_SymStd_UTILS_matvec(k, sod.data(), y.data(),
		                         xd.data(), so.len(0),
		                         so.len(1), kf, ifd, nstencil);
	}
}

}}}
#endif
