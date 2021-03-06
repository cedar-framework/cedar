#ifndef CEDAR_2D_KERNEL_RESIDUAL_H
#define CEDAR_2D_KERNEL_RESIDUAL_H

#include <cedar/kernels/residual.h>
#include <cedar/2d/types.h>

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_residual(int*,real_t*,real_t*,real_t*,real_t*,len_t*,len_t*,
	                          int*,int*,int*,int*,int*,int*,int*);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr2 {

class residual_f90 : public kernels::residual<stypes>
{
	void run(const stencil_op<five_pt> & so,
	         const grid_func & x,
	         const grid_func & b,
	         grid_func & r) override
	{
		this->run_impl(so, x, b, r);
	}
	void run(const stencil_op<nine_pt> & so,
	         const grid_func & x,
	         const grid_func & b,
	         grid_func & r) override
	{
		this->run_impl(so, x, b, r);
	}

	template<class sten>
	void run_impl(const stencil_op<sten> & so,
	              const grid_func & x,
	              const grid_func & b,
	              grid_func & r)
	{
		int k = 0;
		int kf = 0;
		int ifd;
		int ibc;
		int nstncl = stencil_ndirs<sten>::value;
		if (std::is_same<sten, five_pt>::value)
			ifd = 1;
		else
			ifd = 0;
		int irelax = 0;
		int irelax_sym = 0;
		int updown = 0;
		len_t ii = r.len(0);
		len_t jj = r.len(1);

		auto & Ad = const_cast<stencil_op<sten>&>(so);
		grid_func &xd = const_cast<grid_func&>(x);
		grid_func &bd = const_cast<grid_func&>(b);

		BMG_get_bc(params->per_mask(), &ibc);
		BMG2_SymStd_residual(&k, Ad.data(), bd.data(), xd.data(), r.data(), &ii, &jj,
							 &kf, &ifd, &nstncl, &ibc, &irelax, &irelax_sym, &updown);
	}

};

}}

#endif
