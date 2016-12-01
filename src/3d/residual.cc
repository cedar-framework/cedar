#include <boxmg/types.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/grid_func.h>

#include <boxmg/3d/residual.h>

extern "C" {
	using namespace boxmg;
	void BMG3_SymStd_residual(int kg, int NOG, int ifd,
	                          real_t *q, real_t *qf, real_t *so, real_t *RES,
	                          len_t ii, len_t jj, len_t kk,
	                          int NStncl);
}

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	void residual(const stencil_op & A, const grid_func & x,
	              const grid_func & b, grid_func & r)
	{
		int kg, nog, ifd, nstencil;

		auto &Ad = const_cast<stencil_op&>(A);
		auto &xd = const_cast<grid_func&>(x);
		auto &bd = const_cast<grid_func&>(b);

		if (Ad.stencil().five_pt()) {
			ifd = 1;
			nstencil = 4;
			nog = 1;
			kg = 1;
		} else {
			ifd = 0;
			nstencil = 14;
			nog = 1;
			kg = 1;
		}

		BMG3_SymStd_residual(kg, nog, ifd, xd.data(), bd.data(), Ad.data(), r.data(),
		                     r.len(0), r.len(1), r.len(2), nstencil);
	}
}

}}}
