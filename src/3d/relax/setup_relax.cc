#include <boxmg/2d/ftn/BMG_parameters_c.h>
#include <boxmg/3d/relax/setup_relax.h>

extern "C" {
	using namespace boxmg;
	void BMG3_SymStd_SETUP_recip(real_t *so, real_t *sor,
	                             len_t nx, len_t ny, len_t nz,
	                             int nstencl, int nsorv);
}


namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	using namespace boxmg::bmg3;
	void setup_rbgs_point(const stencil_op &so,
	                      relax_stencil &sor)
	{
		int nsorv, nstencil;
		const grid_stencil &sten = so.stencil();
		stencil_op &sod = const_cast<stencil_op&>(so);

		if (sten.five_pt()) {
			nstencil = 4;
		} else {
			nstencil = 14;
		}

		nsorv = 2;

		BMG3_SymStd_SETUP_recip(sod.data(),
		                        sor.data(),
		                        sten.len(0), sten.len(1), sten.len(2),
		                        nstencil, nsorv);
	}


	void setup_relax_xy(const stencil_op &so,
	                    std::vector<slv2_ptr> & planes)
	{
		auto & sten = so.stencil();
		for (auto k : sten.range(2)) {
			auto so2 = bmg2d::stencil_op(sten.shape(0), sten.shape(1));
			auto & sten2 = so2.stencil();
			sten2.five_pt() = sten.five_pt();
			if (sten.five_pt()) {
				for (auto j : sten.range(1)) {
					for (auto i : sten.range(0)) {
						sten2(i,j,bmg2d::dir::C) = sten(i,j,k,dir::P);
						sten2(i,j,bmg2d::dir::W) = sten(i,j,k,dir::PW);
						sten2(i,j,bmg2d::dir::S) = sten(i,j,k,dir::PW);
					}
				}
			} else {
				for (auto j : sten.range(1)) {
					for (auto i : sten.range(0)) {
						sten2(i,j,bmg2d::dir::C) = sten(i,j,k,dir::P);
						sten2(i,j,bmg2d::dir::W) = sten(i,j,k,dir::PW);
						sten2(i,j,bmg2d::dir::S) = sten(i,j,k,dir::PW);
						sten2(i,j,bmg2d::dir::SW) = sten(i,j,k,dir::PSW);
						sten2(i-1,j,bmg2d::dir::SE) = sten(i,j,k,dir::PNW);
					}
				}
			}

			planes.emplace_back(std::make_unique<::boxmg::bmg2d::solver>(std::move(so2)));
			planes.back()->level(-1).x = bmg2d::grid_func(sten.shape(0), sten.shape(1));
			planes.back()->level(-1).b = bmg2d::grid_func(sten.shape(0), sten.shape(1));
		}
	}
}

}}}
