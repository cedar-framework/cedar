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

	void residual(const kernel_params & params,
	              const stencil_op & A, const grid_func & x,
				  const grid_func & b, grid_func &r)
	{
		using namespace cedar::cdr2;

		const grid_stencil &so = A.stencil();

		for (auto j: r.range(1)) {
			for (auto i: r.range(0)) {
				r(i,j) = (b(i,j) +
				          (so(i,j,dir::W)  * x(i-1, j  ) +
				           so(i,j,dir::E)  * x(i+1, j  ) +
				           so(i,j,dir::S)  * x(i  , j-1) +
				           so(i,j,dir::N)  * x(i  , j+1) +
						   so(i,j,dir::SW) * x(i-1, j-1) +
				           so(i,j,dir::SE) * x(i+1, j-1) +
				           so(i,j,dir::NW) * x(i-1, j+1) +
				           so(i,j,dir::NE) * x(i+1, j+1) -
						   so(i,j,dir::C)  * x(i  , j)));
			}
		}
	}

	void residual_fortran(const kernel_params & params,
	                      const stencil_op &A, const grid_func &x,
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

		stencil_op &Ad = const_cast<stencil_op&>(A);
		grid_func &xd = const_cast<grid_func&>(x);
		grid_func &bd = const_cast<grid_func&>(b);
		using namespace cedar::cdr2;

		if (Ad.stencil().five_pt()) {
			ifd = 1;
			nstncl = 3;
		} else {
			ifd = 0;
			nstncl = 5;
		}

		BMG_get_bc(params.per_mask(), &ibc);

		BMG2_SymStd_residual(&k, Ad.data(), bd.data(), xd.data(), r.data(), &ii, &jj,
							 &kf, &ifd, &nstncl, &ibc, &irelax, &irelax_sym, &updown);

	}
}


}}}
