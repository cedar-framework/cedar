#include <boxmg/3d/relax/relax_planes.h>



namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{

	static void copy_rhs_xy(const stencil_op & so,
	                        const grid_func & x,
	                        const grid_func & b,
	                        bmg2d::grid_func & b2,
	                        int ipl)
	{
		auto & o = so.stencil();

		if (o.five_pt()) {
			for (auto j : b.range(1)) {
				for (auto i : b.range(0)) {
					b2(i,j) = b(i,j,ipl)
						+ o(i,j,ipl,dir::B) * x(i,j,ipl-1)
						+ o(i,j,ipl+1,dir::B) * x(i,j,ipl+1);
				}
			}
		} else {
			for (auto j : b.range(1)) {
				for (auto i : b.range(0)) {
					b2(i,j) = b(i,j,ipl)
						+ o(i,j,ipl,dir::B)*x(i,j,ipl-1)
						+ o(i,j,ipl,dir::BW)*x(i-1,j,ipl-1)
						+ o(i,j+1,ipl,dir::BNW)*x(i-1,j+1,ipl-1)
						+ o(i,j+1,ipl,dir::BN)*x(i,j+1,ipl-1)
						+ o(i+1,j+1,ipl,dir::BNE)*x(i+1,j+1,ipl-1)
						+ o(i+1,j,ipl,dir::BE)*x(i+1,j,ipl-1)
						+ o(i+1,j,ipl,dir::BSE)*x(i+1,j-1,ipl-1)
						+ o(i,j,ipl,dir::BS)*x(i,j-1,ipl-1)
						+ o(i,j,ipl,dir::BSW)*x(i-1,j-1,ipl-1)
						+ o(i,j,ipl+1,dir::BE)*x(i-1,j,ipl+1)
						+ o(i,j+1,ipl+1,dir::BSE)*x(i-1,j+1,ipl+1)
						+ o(i,j+1,ipl+1,dir::BS)*x(i,j+1,ipl+1)
						+ o(i+1,j+1,ipl+1,dir::BSW)*x(i+1,j+1,ipl+1)
						+ o(i+1,j,ipl+1,dir::BW)*x(i+1,j,ipl+1)
						+ o(i,j,ipl+1,dir::B)*x(i,j,ipl+1)
						+ o(i+1,j,ipl+1,dir::BNW)*x(i+1,j-1,ipl+1)
						+ o(i,j,ipl+1,dir::BN)*x(i,j-1,ipl+1)
						+ o(i,j,ipl+1,dir::BNE)*x(i-1,j-1,ipl+1);
				}
			}
		}

	}


	static void copy_rhs_xz(const stencil_op & so,
	                        const grid_func & x,
	                        const grid_func & b,
	                        bmg2d::grid_func & b2,
	                        int ipl)
	{
		auto & o = so.stencil();

		if (o.five_pt()) {
			for (auto k : b.range(2)) {
				for (auto i : b.range(0)) {
					b2(i,k) = b(i,ipl,k)
						+ o(i,ipl,k,dir::PS) * x(i,ipl-1,k)
						+ o(i,ipl+1,k,dir::PS) * x(i,ipl+1,k);
				}
			}
		} else {
			for (auto k : b.range(2)) {
				for (auto i : b.range(0)) {
					b2(i,k) = b(i,ipl,k)
						+ o(i,ipl+1,k,dir::PNW) * x(i-1,ipl+1,k)
						+ o(i,ipl+1,k,dir::PS) * x(i,ipl+1,k)
						+ o(i+1,ipl+1,k,dir::PSW) * x(i+1,ipl+1,k)
						+ o(i,ipl+1,k,dir::BNW) * x(i-1,ipl+1,k-1)
						+ o(i,ipl+1,k,dir::BN) * x(i,ipl+1,k-1)
						+ o(i+1,ipl+1,k,dir::BNE) * x(i+1,ipl+1,k-1)
						+ o(i,ipl+1,k+1,dir::BSE) * x(i-1,ipl+1,k+1)
						+ o(i,ipl+1,k+1,dir::BS) * x(i,ipl+1,k+1)
						+ o(i+1,ipl+1,k+1,dir::BSW) * x(i+1,ipl+1,k+1)
						+ o(i,ipl,k,dir::PSW) * x(i-1,ipl-1,k)
						+ o(i,ipl,k,dir::PS) * x(i,ipl-1,k)
						+ o(i+1,ipl,k,dir::PNW) * x(i+1,ipl-1,k)
						+ o(i,ipl,k,dir::BSW) * x(i-1,ipl-1,k-1)
						+ o(i,ipl,k,dir::BS) * x(i,ipl-1,k-1)
						+ o(i+1,ipl,k,dir::BSE) * x(i+1,ipl-1,k-1)
						+ o(i,ipl,k+1,dir::BNE) * x(i-1,ipl-1,k+1)
						+ o(i,ipl,k+1,dir::BN) * x(i,ipl-1,k+1)
						+ o(i+1,ipl,k+1,dir::BNW) * x(i+1,ipl-1,k+1);
				}
			}
		}
	}


	static void copy_rhs_yz(const stencil_op & so,
	                        const grid_func & x,
	                        const grid_func & b,
	                        bmg2d::grid_func & b2,
	                        int ipl)
	{
		auto & o = so.stencil();

		if (o.five_pt()) {
			for (auto k : b.range(2)) {
				for (auto j : b.range(1)) {
					b2(j,k) = b(ipl,j,k)
						+ o(ipl,j,k,dir::PW) * x(ipl-1,j,k)
						+ o(ipl+1,j,k,dir::PW) * x(ipl+1,j,k);
				}
			}
		} else {
			for (auto k : b.range(2)) {
				for (auto j : b.range(1)) {
					b2(j,k) = b(ipl,j,k)
						+ o(ipl,j+1,k,dir::PNW) * x(ipl-1,j+1,k)
						+ o(ipl,j,k,dir::PW) * x(ipl-1,j,k)
						+ o(ipl,j,k,dir::PSW) * x(ipl-1,j-1,k)
						+ o(ipl,j+1,k,dir::BNW) * x(ipl-1,j+1,k-1)
						+ o(ipl,j,k,dir::BW) * x(ipl-1,j,k-1)
						+ o(ipl,j,k,dir::BSW) * x(ipl-1,j-1,k-1)
						+ o(ipl,j+1,k+1,dir::BSE) * x(ipl-1,j+1,k+1)
						+ o(ipl,j,k+1,dir::BE) * x(ipl-1,j,k+1)
						+ o(ipl,j,k+1,dir::BNE) * x(ipl-1,j-1,k+1)
						+ o(ipl+1,j+1,k,dir::PSW) * x(ipl+1,j+1,k)
						+ o(ipl+1,j,k,dir::PW) * x(ipl+1,j,k)
						+ o(ipl+1,j,k,dir::PNW) * x(ipl+1,j-1,k)
						+ o(ipl+1,j+1,k,dir::BNE) * x(ipl+1,j+1,k-1)
						+ o(ipl+1,j,k,dir::BE) * x(ipl+1,j,k-1)
						+ o(ipl+1,j,k,dir::BSE) * x(ipl+1,j-1,k-1)
						+ o(ipl+1,j+1,k+1,dir::BSW) * x(ipl+1,j+1,k+1)
						+ o(ipl+1,j,k+1,dir::BW) * x(ipl+1,j,k+1)
						+ o(ipl+1,j,k+1,dir::BNW) * x(ipl+1,j-1,k+1);
				}
			}
		}
	}


	void relax_xy(const stencil_op & so,
	              grid_func & x,
	              const grid_func & b,
	              cycle::Dir cycle_dir,
	              std::vector<slv2_ptr> & planes)
	{

		auto copy32 = [](grid_func & x, ::boxmg::bmg2d::grid_func & x2, int ipl) {
			for (auto j : x.grange(1)) {
				for (auto i : x.grange(0)) {
					x2(i,j) = x(i,j,ipl);
				}
			}
		};


		auto copy23 = [](::boxmg::bmg2d::grid_func & x2, grid_func & x, int ipl) {
			for (auto j : x.grange(1)) {
				for (auto i : x.grange(0)) {
					x(i,j,ipl) = x2(i,j);
				}
			}
		};


		// indices for symmetric red-black relaxation
		int lstart, lend, lstride;

		if (cycle_dir == cycle::Dir::DOWN) {
			lstart = 1;
			lend = 3;
			lstride = 1;
		} else {
			lstart = 2;
			lend = 0;
			lstride = -1;
		}

		auto nz = x.shape(2);

		// red-black relaxation
		for (auto ipl_beg : range<int>(lstart, lend, lstride)) {
			for (auto ipl : range<int>(ipl_beg, nz+1, 2)) {
				auto & x2 = planes[ipl-1]->level(-1).x;
				auto & b2 = planes[ipl-1]->level(-1).b;

				copy32(x, x2, ipl);
				copy_rhs_xy(so, x, b, b2, ipl);

				log::set_header_msg(" (plane-xy " + std::to_string(ipl) + ")");
				planes[ipl-1]->solve(b2, x2);
				log::set_header_msg("");

				copy23(x2, x, ipl);
			}
		}
	}


	void relax_xz(const stencil_op & so,
	              grid_func & x,
	              const grid_func & b,
	              cycle::Dir cycle_dir,
	              std::vector<slv2_ptr> & planes)
	{

		auto copy32 = [](grid_func & x, ::boxmg::bmg2d::grid_func & x2, int ipl) {
			for (auto k : x.grange(2)) {
				for (auto i : x.grange(0)) {
					x2(i,k) = x(i,ipl,k);
				}
			}
		};


		auto copy23 = [](::boxmg::bmg2d::grid_func & x2, grid_func & x, int ipl) {
			for (auto k : x.grange(2)) {
				for (auto i : x.grange(0)) {
					x(i,ipl,k) = x2(i,k);
				}
			}
		};


		// indices for symmetric red-black relaxation
		int lstart, lend, lstride;

		if (cycle_dir == cycle::Dir::DOWN) {
			lstart = 1;
			lend = 3;
			lstride = 1;
		} else {
			lstart = 2;
			lend = 0;
			lstride = -1;
		}

		auto ny = x.shape(1);

		// red-black relaxation
		for (auto ipl_beg : range<int>(lstart, lend, lstride)) {
			for (auto ipl : range<int>(ipl_beg, ny+1, 2)) {
				auto & x2 = planes[ipl-1]->level(-1).x;
				auto & b2 = planes[ipl-1]->level(-1).b;

				copy32(x, x2, ipl);
				copy_rhs_xz(so, x, b, b2, ipl);

				log::set_header_msg(" (plane-xz " + std::to_string(ipl) + ")");
				planes[ipl-1]->solve(b2, x2);
				log::set_header_msg("");

				copy23(x2, x, ipl);
			}
		}
	}


	void relax_yz(const stencil_op & so,
	              grid_func & x,
	              const grid_func & b,
	              cycle::Dir cycle_dir,
	              std::vector<slv2_ptr> & planes)
	{

		auto copy32 = [](grid_func & x, ::boxmg::bmg2d::grid_func & x2, int ipl) {
			for (auto k : x.grange(2)) {
				for (auto j : x.grange(1)) {
					x2(j,k) = x(ipl,j,k);
				}
			}
		};


		auto copy23 = [](::boxmg::bmg2d::grid_func & x2, grid_func & x, int ipl) {
			for (auto k : x.grange(2)) {
				for (auto j : x.grange(1)) {
					x(ipl,j,k) = x2(j,k);
				}
			}
		};


		// indices for symmetric red-black relaxation
		int lstart, lend, lstride;

		if (cycle_dir == cycle::Dir::DOWN) {
			lstart = 1;
			lend = 3;
			lstride = 1;
		} else {
			lstart = 2;
			lend = 0;
			lstride = -1;
		}

		auto nx = x.shape(0);

		// red-black relaxation
		for (auto ipl_beg : range<int>(lstart, lend, lstride)) {
			for (auto ipl : range<int>(ipl_beg, nx+1, 2)) {
				auto & x2 = planes[ipl-1]->level(-1).x;
				auto & b2 = planes[ipl-1]->level(-1).b;

				copy32(x, x2, ipl);
				copy_rhs_yz(so, x, b, b2, ipl);

				log::set_header_msg(" (plane-yz " + std::to_string(ipl) + ")");
				planes[ipl-1]->solve(b2, x2);
				log::set_header_msg("");

				copy23(x2, x, ipl);
			}
		}
	}
}

}}}
