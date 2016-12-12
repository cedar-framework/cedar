#include <boxmg/2d/kernel/factory.h>
#include <boxmg/2d/gallery.h>

namespace boxmg { namespace bmg2d {

using namespace boxmg;

stencil_op create_poisson(len_t nx, len_t ny)
{
	config::reader conf("");
	auto kreg = kernel::factory::from_config(conf);

	auto so = stencil_op(nx, ny);
	so.set_registry(kreg);
	auto & sten = so.stencil();
	sten.five_pt() = true;

	sten.set(0);

	real_t hx = 1.0/(sten.len(0)-1);
	real_t hy = 1.0/(sten.len(1)-1);
	real_t xh = hy/hx;
	real_t yh = hx/hy;
	len_t i1 = sten.shape(0)+1;
	len_t j1 = sten.shape(1)+1;

	for (auto j : range<len_t>(2, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			sten(i,j,dir::S) = 1.0 * yh;
		}
	}

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(2, i1)) {
			sten(i,j,dir::W) = 1.0 * xh;
		}
	}

	for (auto j : sten.range(1)) {
		for (auto i : sten.range(0)) {
			sten(i,j,dir::C) = 2*xh + 2*yh;
		}
	}

	return so;
}


stencil_op create_diag_diffusion(len_t nx, len_t ny, real_t dx, real_t dy)
{
	config::reader conf("");
	auto kreg = kernel::factory::from_config(conf);

	auto so = stencil_op(nx, ny);
	so.set_registry(kreg);
	auto & sten = so.stencil();
	sten.five_pt() = true;

	sten.set(0);

	real_t hx = 1.0/(sten.len(0)-1);
	real_t hy = 1.0/(sten.len(1)-1);
	real_t xh = hy/hx;
	real_t yh = hx/hy;

	len_t i1 = sten.shape(0)+1;
	len_t j1 = sten.shape(1)+1;

	for (auto j : range<len_t>(2, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			sten(i,j,dir::S) = dy*yh;
		}
	}

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(2, i1)) {
			sten(i,j,dir::W) = dx*xh;
		}
	}

	for (auto j : sten.range(1)) {
		for (auto i : sten.range(0)) {
			sten(i,j,dir::C) = 2*dx*xh + 2*dy*yh;
		}
	}

	return so;
}


stencil_op create_fe(len_t nx, len_t ny)
{
	config::reader conf("");
	auto kreg = kernel::factory::from_config(conf);

	auto so = stencil_op(nx, ny);
	so.set_registry(kreg);
	auto & sten = so.stencil();
	sten.five_pt() = false;

	sten.set(0);

	len_t i1 = sten.shape(0)+1;
	len_t j1 = sten.shape(1)+1;

	for (auto j : range<len_t>(2, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			sten(i,j,dir::S) = 1.0;
		}
	}

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(2, i1)) {
			sten(i,j,dir::W) = 1.0;
		}
	}

	for (auto j : range<len_t>(2, j1)) {
		for (auto i : range<len_t>(2, i1)) {
			sten(i,j,dir::SW) = 1.0;
			sten(i-1,j,dir::SE) = 1.0;
		}
	}

	for (auto j: sten.range(1)) {
		for (auto i : sten.range(0)) {
			sten(i,j,dir::C) = 8.0;
		}
	}

	return so;
}

}}
