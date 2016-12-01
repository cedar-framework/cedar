#include "boxmg/2d/residual.h"

	extern "C" {
		using namespace boxmg;
		void BMG2_SymStd_residual(int*,real_t*,real_t*,real_t*,real_t*,len_t*,len_t*,
								  int*,int*,int*,int*,int*,int*,int*);
	}

namespace boxmg { namespace bmg2d { namespace kernel {


namespace impls
{
	using namespace boxmg::bmg2d;

	void residual(const stencil_op & A, const grid_func & x,
				  const grid_func & b, grid_func &r)
	{
		using namespace boxmg::bmg2d;

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

	void residual_fortran(const stencil_op &A, const grid_func &x,
						  const grid_func &b, grid_func &r)
	{
		int k = 0;
		int kf = 0;
		int ifd = 0;
		int nstncl = 5;
		int ibc = 0;
		int irelax = 0;
		int irelax_sym = 0;
		int updown = 0;
		len_t ii = r.len(0);
		len_t jj = r.len(1);

		stencil_op &Ad = const_cast<stencil_op&>(A);
		grid_func &xd = const_cast<grid_func&>(x);
		grid_func &bd = const_cast<grid_func&>(b);
		using namespace boxmg::bmg2d;

		BMG2_SymStd_residual(&k, Ad.data(), bd.data(), xd.data(), r.data(), &ii, &jj,
							 &kf, &ifd, &nstncl, &ibc, &irelax, &irelax_sym, &updown);

	}
}


}}}
