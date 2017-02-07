#include "boxmg/2d/ftn/BMG_parameters_c.h"
#include "boxmg/2d/relax/relax.h"

extern "C" {
	using namespace boxmg;
	void BMG2_SymStd_relax_GS(int, real_t*, real_t*, real_t*, real_t*, len_t, len_t,
	                          int, int, int, int, int, int, int);
	void BMG2_SymStd_relax_lines_x(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                               real_t *B, len_t II, len_t JJ, int kf, int ifd,
	                               int nstencil, int irelax_sym, int updown, int jpn);
	void BMG2_SymStd_relax_lines_y(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                               real_t *B, len_t II, len_t JJ, int kf, int ifd,
	                               int nstencil, int irelax_sym, int updown, int jpn);
}

namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	void relax_rbgs_point(const stencil_op & so,
	                      grid_func & x,
	                      const grid_func & b,
	                      const relax_stencil & sor,
	                      cycle::Dir cycle_dir)
	{
		using namespace boxmg::bmg2d;
		int k, kf, ifd;
		int updown, nsorv, ibc, nstencil;

		const grid_stencil &so_sten = so.stencil();
		stencil_op & sod = const_cast<stencil_op&>(so);
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		grid_func & bd = const_cast<grid_func&>(b);

		k = kf = 1;
		if (so_sten.five_pt()) {
			ifd = 1;
			nstencil = 3;
		} else {
			ifd = 0;
			nstencil = 5;
		}

		nsorv = 2;

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		ibc = BMG_BCs_definite;

		BMG2_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                     so_sten.len(0), so_sten.len(1), kf, ifd, nstencil, nsorv,
		                     BMG_RELAX_SYM, updown, ibc);
	}



	void relax_lines_x(const stencil_op & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func & res,
	                   cycle::Dir cycle_dir)
	{
		using namespace boxmg::bmg2d;
		int k, kf, ifd;
		int updown, nsorv, ibc, nstencil;

		const grid_stencil &so_sten = so.stencil();
		stencil_op & sod = const_cast<stencil_op&>(so);
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		grid_func & bd = const_cast<grid_func&>(b);

		k = kf = 1;
		if (so_sten.five_pt()) {
			ifd = 1;
			nstencil = 3;
		} else {
			ifd = 0;
			nstencil = 5;
		}

		nsorv = 2;

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		ibc = BMG_BCs_definite;

		BMG2_SymStd_relax_lines_x(k, sod.data(), bd.data(), x.data(), sord.data(), res.data(),
		                          so_sten.len(0), so_sten.len(1), kf, ifd, nstencil,
		                          BMG_RELAX_SYM, updown, ibc);
	}


	void relax_lines_y(const stencil_op & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func & res,
	                   cycle::Dir cycle_dir)
	{
		using namespace boxmg::bmg2d;
		int k, kf, ifd;
		int updown, nsorv, ibc, nstencil;

		const grid_stencil &so_sten = so.stencil();
		stencil_op & sod = const_cast<stencil_op&>(so);
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		grid_func & bd = const_cast<grid_func&>(b);

		k = kf = 1;
		if (so_sten.five_pt()) {
			ifd = 1;
			nstencil = 3;
		} else {
			ifd = 0;
			nstencil = 5;
		}

		nsorv = 2;

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		ibc = BMG_BCs_definite;

		BMG2_SymStd_relax_lines_y(k, sod.data(), bd.data(), x.data(), sord.data(), res.data(),
		                          so_sten.len(0), so_sten.len(1), kf, ifd, nstencil,
		                          BMG_RELAX_SYM, updown, ibc);
	}
}

}}}
