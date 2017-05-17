#include <cedar/2d/gallery.h>

namespace cedar { namespace cdr2 { namespace gallery {

using namespace cedar;

stencil_op<five_pt> poisson(len_t nx, len_t ny)
{
	auto so = stencil_op<five_pt>(nx, ny);

	so.set(0);

	real_t hx = 1.0/(so.len(0)-1);
	real_t hy = 1.0/(so.len(1)-1);
	real_t xh = hy/hx;
	real_t yh = hx/hy;
	len_t i1 = so.shape(0)+1;
	len_t j1 = so.shape(1)+1;

	for (auto j : range<len_t>(2, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			so(i,j,five_pt::s) = 1.0 * yh;
		}
	}

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(2, i1)) {
			so(i,j,five_pt::w) = 1.0 * xh;
		}
	}

	for (auto j : so.range(1)) {
		for (auto i : so.range(0)) {
			so(i,j,five_pt::c) = 2*xh + 2*yh;
		}
	}

	return so;
}


stencil_op<five_pt> diag_diffusion(len_t nx, len_t ny, real_t dx, real_t dy)
{
	auto so = stencil_op<five_pt>(nx, ny);

	so.set(0);

	real_t hx = 1.0/(so.len(0)-1);
	real_t hy = 1.0/(so.len(1)-1);
	real_t xh = hy/hx;
	real_t yh = hx/hy;

	len_t i1 = so.shape(0)+1;
	len_t j1 = so.shape(1)+1;

	for (auto j : range<len_t>(2, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			so(i,j,five_pt::s) = dy*yh;
		}
	}

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(2, i1)) {
			so(i,j,five_pt::w) = dx*xh;
		}
	}

	for (auto j : so.range(1)) {
		for (auto i : so.range(0)) {
			so(i,j,five_pt::c) = 2*dx*xh + 2*dy*yh;
		}
	}

	return so;
}


stencil_op<nine_pt> fe(len_t nx, len_t ny)
{
	auto so = stencil_op<nine_pt>(nx, ny);

	so.set(0);

	len_t i1 = so.shape(0)+1;
	len_t j1 = so.shape(1)+1;

	for (auto j : range<len_t>(2, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			so(i,j,nine_pt::s) = 1.0;
		}
	}

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(2, i1)) {
			so(i,j,nine_pt::w) = 1.0;
		}
	}

	for (auto j : range<len_t>(2, j1)) {
		for (auto i : range<len_t>(2, i1)) {
			so(i,j,nine_pt::sw) = 1.0;
			so(i,j,nine_pt::nw) = 1.0;
		}
	}

	for (auto j: so.range(1)) {
		for (auto i : so.range(0)) {
			so(i,j,nine_pt::c) = 8.0;
		}
	}

	return so;
}

}}}
