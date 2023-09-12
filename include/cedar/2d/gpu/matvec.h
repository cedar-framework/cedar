#ifndef CEDAR_2D_GPU_MATVEC_H
#define CEDAR_2D_GPU_MATVEC_H

#include <cedar/2d/gpu/types.h>
#include <cedar/kernels/matvec.h>

using real_t = cedar::real_t;
using len_t = cedar::len_t;
#include <src/2d/ftn/mpi/BMG2_SymStd_UTILS_matvec.f90.hpp>

namespace cedar { namespace cdr2 { namespace gpu { namespace mpi {

class matvec_f90 : public kernels::matvec<stypes>
{
	void run(const stencil_op<five_pt> & so,
	         const grid_func & x,
	         grid_func & y) override
	{ this->run_impl(so, x, y); }
	void run(const stencil_op<nine_pt> & so,
	         const grid_func & x,
	         grid_func & y) override
	{ this->run_impl(so, x, y); }

	template<class sten>
	void run_impl(const stencil_op<sten> & so,
	              const grid_func & x,
	              grid_func & y)
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

		BMG2_SymStd_UTILS_matvec<ftl::device::GPU>(
                    k, sod, y, xd, so.len(0), so.len(1), kf, ifd, nstencil);
	}
};

}}}}

#endif
