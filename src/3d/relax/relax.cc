#include "boxmg/2d/ftn/BMG_parameters_c.h"

#include <boxmg/3d/relax/relax.h>

extern "C" {
	using namespace boxmg;
	void BMG3_SymStd_relax_GS(int kg, real_t *so, real_t *qf, real_t *q, real_t *sor,
	                          len_t ii, len_t jj, len_t kk, int ifd, int nstncl, int nsorv,
	                          int irelax_sym, int updown, int jpn);
}


namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void relax_rbgs_point(const stencil_op & so,
	                      grid_func &x,
	                      const grid_func &b,
	                      const relax_stencil & sor,
	                      cycle::Dir cycle_dir)
	{
		using namespace boxmg::bmg3;
		int k, ifd, nstencil, nsorv, jpn, updown;

		const grid_stencil &so_sten = so.stencil();
		auto & sod = const_cast<stencil_op&>(so);
		auto & sord = const_cast<relax_stencil&>(sor);
		auto & bd = const_cast<grid_func&>(b);

		nsorv = 2;
		// give these dummy values, ifd will pass the relevant information.
		k = 1;

		if (so_sten.five_pt()) {
			ifd = 1;
			nstencil = 4;
		} else {
			ifd = 0;
			nstencil = 14;
		}

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		jpn = BMG_BCs_definite;

		BMG3_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                     so_sten.len(0), so_sten.len(1), so_sten.len(2), ifd, nstencil, nsorv,
		                     BMG_RELAX_SYM, updown, jpn);
	}

	void relax_xy(const stencil_op & so,
	              grid_func & x,
	              const grid_func & b,
	              cycle::Dir cycle_dir,
	              std::vector<slv2_ptr> & planes)
	{
	}
}

}}}
