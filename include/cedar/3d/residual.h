#ifndef CEDAR_3D_KERNEL_RESIDUAL_H
#define CEDAR_3D_KERNEL_RESIDUAL_H

#include <type_traits>
#include <cedar/3d/types.h>
#include <cedar/kernels/residual.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_residual(int kg, int NOG, int ifd,
	                          real_t *q, real_t *qf, real_t *so, real_t *RES,
	                          len_t ii, len_t jj, len_t kk,
	                          int NStncl);
}


namespace cedar { namespace cdr3 {
class residual_f90 : public kernels::residual<stypes>
{
	void run(const stencil_op<seven_pt> & so,
	         const grid_func & x,
	         const grid_func & b,
	         grid_func & r) override
	{
		this->run_impl(so, x, b, r);
	}
	void run(const stencil_op<xxvii_pt> & so,
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
		int kg, nog, ifd, nstencil;

		auto &Ad = const_cast<stencil_op<sten>&>(so);
		auto &xd = const_cast<grid_func&>(x);
		auto &bd = const_cast<grid_func&>(b);

		nstencil = stencil_ndirs<sten>::value;
		if (std::is_same<sten,seven_pt>::value) {
			ifd = 1;
			nog = 1;
			kg = 1;
		} else {
			ifd = 0;
			nog = 1;
			kg = 1;
		}

		BMG3_SymStd_residual(kg, nog, ifd, xd.data(), bd.data(), Ad.data(), r.data(),
		                     r.len(0), r.len(1), r.len(2), nstencil);
	}
};

}}

#endif
