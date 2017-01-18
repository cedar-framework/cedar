#include <boxmg/3d/relax/setup_planes.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	void setup_relax_xy(const stencil_op &so,
	                    std::vector<slv2_ptr> & planes)
	{
		auto & sten = so.stencil();
		for (auto k : sten.range(2)) {
			auto so2 = bmg2d::stencil_op(sten.shape(0), sten.shape(1));
			auto & sten2 = so2.stencil();
			sten2.five_pt() = sten.five_pt();
			if (sten.five_pt()) {
				for (auto j : sten.grange(1)) {
					for (auto i : sten.grange(0)) {
						sten2(i,j,bmg2d::dir::C) = sten(i,j,k,dir::P);
						sten2(i,j,bmg2d::dir::W) = sten(i,j,k,dir::PW);
						sten2(i,j,bmg2d::dir::S) = sten(i,j,k,dir::PS);
					}
				}
			} else {
				for (auto j : sten.grange(1)) {
					for (auto i : sten.grange(0)) {
						sten2(i,j,bmg2d::dir::C) = sten(i,j,k,dir::P);
						sten2(i,j,bmg2d::dir::W) = sten(i,j,k,dir::PW);
						sten2(i,j,bmg2d::dir::S) = sten(i,j,k,dir::PS);
						sten2(i,j,bmg2d::dir::SW) = sten(i,j,k,dir::PSW);
						sten2(i-1,j,bmg2d::dir::SE) = sten(i,j,k,dir::PNW);
					}
				}
			}

			auto conf2 = config::reader("config.json");
			conf2.set("solver.relaxation", "line-xy");

			planes.emplace_back(std::make_unique<::boxmg::bmg2d::solver>(std::move(so2), std::move(conf2)));
			planes.back()->level(-1).x = bmg2d::grid_func(sten.shape(0), sten.shape(1));
			planes.back()->level(-1).b = bmg2d::grid_func(sten.shape(0), sten.shape(1));
		}
	}


	void setup_relax_xz(const stencil_op &so,
	                    std::vector<slv2_ptr> & planes)
	{
		auto & o = so.stencil();

		for (auto j : o.range(1)) {
			auto so2 = bmg2d::stencil_op(o.shape(0), o.shape(2));
			auto & o2 = so2.stencil();
			o2.five_pt() = o.five_pt();

			if (o.five_pt()) {
				for (auto k : o.grange(2)) {
					for (auto i : o.grange(0)) {
						o2(i,k,bmg2d::dir::C) = o(i,j,k,dir::P);
						o2(i,k,bmg2d::dir::W) = o(i,j,k,dir::PW);
						o2(i,k,bmg2d::dir::S) = o(i,j,k,dir::B);
					}
				}
			} else {
				for (auto k : o.grange(2)) {
					for (auto i : o.grange(0)) {
						o2(i,k,bmg2d::dir::C) = o(i,j,k,dir::P);
						o2(i,k,bmg2d::dir::W) = o(i,j,k,dir::PW);
						o2(i,k,bmg2d::dir::S) = o(i,j,k,dir::B);
						o2(i,k,bmg2d::dir::SW) = o(i,j,k,dir::BW);
						o2(i-1,k,bmg2d::dir::SE) = o(i,j,k,dir::BE);
					}
				}
			}

			auto conf2 = config::reader("config.json");
			conf2.set("solver.relaxation", "line-xy");

			planes.emplace_back(std::make_unique<::boxmg::bmg2d::solver>(std::move(so2), std::move(conf2)));
			planes.back()->level(-1).x = bmg2d::grid_func(o.shape(0), o.shape(2));
			planes.back()->level(-1).b = bmg2d::grid_func(o.shape(0), o.shape(2));
		}
	}


	void setup_relax_yz(const stencil_op &so,
	                    std::vector<slv2_ptr> & planes)
	{
		auto & o = so.stencil();

		for (auto i : o.range(0)) {
			auto so2 = bmg2d::stencil_op(o.shape(1), o.shape(2));
			auto & o2 = so2.stencil();
			o2.five_pt() = o.five_pt();

			if (o.five_pt()) {
				for (auto k : o.grange(2)) {
					for (auto j : o.grange(1)) {
						o2(j,k,bmg2d::dir::C) = o(i,j,k,dir::P);
						o2(j,k,bmg2d::dir::W) = o(i,j,k,dir::PS);
						o2(j,k,bmg2d::dir::S) = o(i,j,k,dir::B);
					}
				}
			} else {
				for (auto k : o.grange(2)) {
					for (auto j : o.grange(1)) {
						o2(j,k,bmg2d::dir::C) = o(i,j,k,dir::P);
						o2(j,k,bmg2d::dir::W) = o(i,j,k,dir::PS);
						o2(j,k,bmg2d::dir::S) = o(i,j,k,dir::B);
						o2(j,k,bmg2d::dir::SW) = o(i,j,k,dir::BS);
						o2(j-1,k,bmg2d::dir::SE) = o(i,j,k,dir::BN);
					}
				}
			}

			auto conf2 = config::reader("config.json");
			conf2.set("solver.relaxation", "line-xy");

			planes.emplace_back(std::make_unique<::boxmg::bmg2d::solver>(std::move(so2), std::move(conf2)));
			planes.back()->level(-1).x = bmg2d::grid_func(o.shape(1), o.shape(2));
			planes.back()->level(-1).b = bmg2d::grid_func(o.shape(1), o.shape(2));
		}
	}
}

}}}
