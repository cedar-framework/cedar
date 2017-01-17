#include <boxmg/3d/kernel/factory.h>
#include <boxmg/3d/gallery.h>

namespace boxmg { namespace bmg3 { namespace gallery {

using namespace boxmg;

stencil_op poisson(len_t nx, len_t ny, len_t nz)
{
	config::reader conf("");
	auto kreg = kernel::factory::from_config(conf);

	auto so = stencil_op(nx, ny, nz);
	so.set_registry(kreg);
	auto & sten = so.stencil();
	sten.five_pt() = true;

	sten.set(0);

	real_t hx = 1.0/(sten.len(0)-1);
	real_t hy = 1.0/(sten.len(1)-1);
	real_t hz = 1.0/(sten.len(2)-1);
	real_t xh=hy*hz/hx;
	real_t yh=hx*hz/hy;
	real_t zh=hx*hy/hz;
	len_t i1 = sten.shape(0) + 1;
	len_t j1 = sten.shape(1) + 1;
	len_t k1 = sten.shape(2) + 1;

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(2, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				sten(i,j,k,dir::PS) = 1.0*yh;
			}
		}
	}

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(2, i1)) {
				sten(i,j,k,dir::PW) = 1.0*xh;
			}
		}
	}

	for (auto k : range<len_t>(2, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				sten(i,j,k,dir::B) = 1.0*zh;
			}
		}
	}

	for (auto k : sten.range(2)) {
		for (auto j : sten.range(1)) {
			for (auto i : sten.range(0)) {
				sten(i,j,k,dir::P) = 2.0*xh + 2.0*yh + 2.0*zh;
			}
		}
	}

	return so;
}

stencil_op diag_diffusion(len_t nx, len_t ny, len_t nz,
                          real_t dx, real_t dy, real_t dz)
{
	config::reader conf("");
	auto kreg = kernel::factory::from_config(conf);

	auto so = stencil_op(nx, ny, nz);
	so.set_registry(kreg);
	auto & sten = so.stencil();
	sten.five_pt() = true;

	sten.set(0);

	real_t hx = 1.0/(sten.len(0)-1);
	real_t hy = 1.0/(sten.len(1)-1);
	real_t hz = 1.0/(sten.len(2)-1);
	real_t xh=hy*hz/hx;
	real_t yh=hx*hz/hy;
	real_t zh=hx*hy/hz;
	len_t i1 = sten.shape(0) + 1;
	len_t j1 = sten.shape(1) + 1;
	len_t k1 = sten.shape(2) + 1;

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(2, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				sten(i,j,k,dir::PS) = dy*yh;
			}
		}
	}

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(2, i1)) {
				sten(i,j,k,dir::PW) = dx*xh;
			}
		}
	}

	for (auto k : range<len_t>(2, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				sten(i,j,k,dir::B) = dz*zh;
			}
		}
	}

	for (auto k : sten.range(2)) {
		for (auto j : sten.range(1)) {
			for (auto i : sten.range(0)) {
				sten(i,j,k,dir::P) = 2.0*dx*xh + 2.0*dy*yh + 2.0*dz*zh;
			}
		}
	}

	return so;
}

}}}
