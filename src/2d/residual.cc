#include "cedar/2d/residual.h"

	extern "C" {
		using namespace cedar;
		void BMG2_SymStd_residual(int*,real_t*,real_t*,real_t*,real_t*,len_t*,len_t*,
								  int*,int*,int*,int*,int*,int*,int*);
		void BMG_get_bc(int, int*);
	}

namespace cedar { namespace cdr2 { namespace kernel {


namespace impls
{
	using namespace cedar::cdr2;

	template<>
	void residual_fortran(const kernel_params & params,
	                      const stencil_op<five_pt> &A, const grid_func &x,
						  const grid_func &b, grid_func &r)
	{
		int k = 0;
		int kf = 0;
		int ifd = 0;
		int ibc;
		int nstncl = 5;
		int irelax = 0;
		int irelax_sym = 0;
		int updown = 0;
		len_t ii = r.len(0);
		len_t jj = r.len(1);

		auto & Ad = const_cast<stencil_op<five_pt>&>(A);
		grid_func &xd = const_cast<grid_func&>(x);
		grid_func &bd = const_cast<grid_func&>(b);
		using namespace cedar::cdr2;

		ifd = 1;
		nstncl = 3;

		BMG_get_bc(params.per_mask(), &ibc);

		BMG2_SymStd_residual(&k, Ad.data(), bd.data(), xd.data(), r.data(), &ii, &jj,
							 &kf, &ifd, &nstncl, &ibc, &irelax, &irelax_sym, &updown);

	}

	template<>
	void residual_fortran(const kernel_params & params,
	                      const stencil_op<nine_pt> &A, const grid_func &x,
						  const grid_func &b, grid_func &r)
	{
		int k = 0;
		int kf = 0;
		int ifd = 0;
		int ibc;
		int nstncl = 5;
		int irelax = 0;
		int irelax_sym = 0;
		int updown = 0;
		len_t ii = r.len(0);
		len_t jj = r.len(1);

		auto & Ad = const_cast<stencil_op<nine_pt>&>(A);
		grid_func &xd = const_cast<grid_func&>(x);
		grid_func &bd = const_cast<grid_func&>(b);
		using namespace cedar::cdr2;

		ifd = 0;
		nstncl = 5;

		BMG_get_bc(params.per_mask(), &ibc);

		BMG2_SymStd_residual(&k, Ad.data(), bd.data(), xd.data(), r.data(), &ii, &jj,
							 &kf, &ifd, &nstncl, &ibc, &irelax, &irelax_sym, &updown);

	}
}


}}}
