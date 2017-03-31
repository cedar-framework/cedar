#include <cedar/3d/relax/setup_planes.h>

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	static std::shared_ptr<config::reader> plane_config()
	{
		auto conf = std::make_shared<config::reader>("plane.json");
		conf->set("solver.relaxation", "line-xy");
		conf->set("solver.max-iter", 1);

		return conf;
	}


	void setup_relax_xy(const stencil_op &so,
	                    std::vector<slv2_ptr> & planes)
	{
		auto & sten = so.stencil();
		for (auto k : sten.range(2)) {
			auto so2 = cdr2::stencil_op(sten.shape(0), sten.shape(1));
			auto & sten2 = so2.stencil();
			sten2.five_pt() = sten.five_pt();
			if (sten.five_pt()) {
				for (auto j : sten.grange(1)) {
					for (auto i : sten.grange(0)) {
						sten2(i,j,cdr2::dir::C) = sten(i,j,k,dir::P);
						sten2(i,j,cdr2::dir::W) = sten(i,j,k,dir::PW);
						sten2(i,j,cdr2::dir::S) = sten(i,j,k,dir::PS);
					}
				}
			} else {
				for (auto j : sten.grange(1)) {
					for (auto i : sten.grange(0)) {
						sten2(i,j,cdr2::dir::C) = sten(i,j,k,dir::P);
						sten2(i,j,cdr2::dir::W) = sten(i,j,k,dir::PW);
						sten2(i,j,cdr2::dir::S) = sten(i,j,k,dir::PS);
						sten2(i,j,cdr2::dir::SW) = sten(i,j,k,dir::PSW);
						sten2(i-1,j,cdr2::dir::SE) = sten(i,j,k,dir::PNW);
					}
				}
			}

			auto conf2 = plane_config();

			planes.emplace_back(std::make_unique<::cedar::cdr2::solver>(std::move(so2), conf2));
			planes.back()->level(-1).x = cdr2::grid_func(sten.shape(0), sten.shape(1));
			planes.back()->level(-1).b = cdr2::grid_func(sten.shape(0), sten.shape(1));
		}
	}


	void setup_relax_xz(const stencil_op &so,
	                    std::vector<slv2_ptr> & planes)
	{
		auto & o = so.stencil();

		for (auto j : o.range(1)) {
			auto so2 = cdr2::stencil_op(o.shape(0), o.shape(2));
			auto & o2 = so2.stencil();
			o2.five_pt() = o.five_pt();

			if (o.five_pt()) {
				for (auto k : o.grange(2)) {
					for (auto i : o.grange(0)) {
						o2(i,k,cdr2::dir::C) = o(i,j,k,dir::P);
						o2(i,k,cdr2::dir::W) = o(i,j,k,dir::PW);
						o2(i,k,cdr2::dir::S) = o(i,j,k,dir::B);
					}
				}
			} else {
				for (auto k : o.grange(2)) {
					for (auto i : o.grange(0)) {
						o2(i,k,cdr2::dir::C) = o(i,j,k,dir::P);
						o2(i,k,cdr2::dir::W) = o(i,j,k,dir::PW);
						o2(i,k,cdr2::dir::S) = o(i,j,k,dir::B);
						o2(i,k,cdr2::dir::SW) = o(i,j,k,dir::BW);
						o2(i-1,k,cdr2::dir::SE) = o(i,j,k,dir::BE);
					}
				}
			}

			auto conf2 = plane_config();

			planes.emplace_back(std::make_unique<::cedar::cdr2::solver>(std::move(so2), conf2));
			planes.back()->level(-1).x = cdr2::grid_func(o.shape(0), o.shape(2));
			planes.back()->level(-1).b = cdr2::grid_func(o.shape(0), o.shape(2));
		}
	}


	void setup_relax_yz(const stencil_op &so,
	                    std::vector<slv2_ptr> & planes)
	{
		auto & o = so.stencil();

		for (auto i : o.range(0)) {
			auto so2 = cdr2::stencil_op(o.shape(1), o.shape(2));
			auto & o2 = so2.stencil();
			o2.five_pt() = o.five_pt();

			if (o.five_pt()) {
				for (auto k : o.grange(2)) {
					for (auto j : o.grange(1)) {
						o2(j,k,cdr2::dir::C) = o(i,j,k,dir::P);
						o2(j,k,cdr2::dir::W) = o(i,j,k,dir::PS);
						o2(j,k,cdr2::dir::S) = o(i,j,k,dir::B);
					}
				}
			} else {
				for (auto k : o.grange(2)) {
					for (auto j : o.grange(1)) {
						o2(j,k,cdr2::dir::C) = o(i,j,k,dir::P);
						o2(j,k,cdr2::dir::W) = o(i,j,k,dir::PS);
						o2(j,k,cdr2::dir::S) = o(i,j,k,dir::B);
						o2(j,k,cdr2::dir::SW) = o(i,j,k,dir::BS);
						o2(j-1,k,cdr2::dir::SE) = o(i,j,k,dir::BN);
					}
				}
			}

			auto conf2 = plane_config();

			planes.emplace_back(std::make_unique<::cedar::cdr2::solver>(std::move(so2), conf2));
			planes.back()->level(-1).x = cdr2::grid_func(o.shape(1), o.shape(2));
			planes.back()->level(-1).b = cdr2::grid_func(o.shape(1), o.shape(2));
		}
	}
}

}}}
