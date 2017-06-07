#ifndef CEDAR_3D_RELAX_H
#define CEDAR_3D_RELAX_H

#include <type_traits>

#include "cedar/2d/ftn/BMG_parameters_c.h"

#include <cedar/kernel_params.h>
#include <cedar/cycle/types.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/grid_func.h>
#include <cedar/3d/relax_stencil.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_relax_GS(int kg, real_t *so, real_t *qf, real_t *q, real_t *sor,
	                          len_t ii, len_t jj, len_t kk, int ifd, int nstncl, int nsorv,
	                          int irelax_sym, int updown, int jpn);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	template<class sten>
	void relax_rbgs_point(const kernel_params & params,
	                      const stencil_op<sten> & so,
	                      grid_func &x,
	                      const grid_func &b,
	                      const relax_stencil & sor,
	                      cycle::Dir cycle_dir)
	{
		using namespace cedar::cdr3;
		int k, ifd, nstencil, nsorv, jpn, updown;

		auto & sod = const_cast<stencil_op<sten>&>(so);
		auto & sord = const_cast<relax_stencil&>(sor);
		auto & bd = const_cast<grid_func&>(b);

		nsorv = 2;
		// give these dummy values, ifd will pass the relevant information.
		k = 1;

		nstencil = stencil_ndirs<sten>::value;
		if (std::is_same<sten, seven_pt>::value)
			ifd = 1;
		else
			ifd = 0;

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		BMG_get_bc(params.per_mask(), &jpn);

		BMG3_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                     so.len(0), so.len(1), so.len(2), ifd, nstencil, nsorv,
		                     BMG_RELAX_SYM, updown, jpn);
	}
}

}}}

#endif
