#include <cedar/3d/gallery.h>

namespace cedar { namespace cdr3 { namespace gallery {

using namespace cedar;

stencil_op<seven_pt> poisson(len_t nx, len_t ny, len_t nz)
{
	stencil_op<seven_pt> so(nx, ny, nz);

	so.set(0);

	real_t hx = 1.0/(so.len(0)-1);
	real_t hy = 1.0/(so.len(1)-1);
	real_t hz = 1.0/(so.len(2)-1);
	real_t xh=hy*hz/hx;
	real_t yh=hx*hz/hy;
	real_t zh=hx*hy/hz;
	len_t i1 = so.shape(0) + 1;
	len_t j1 = so.shape(1) + 1;
	len_t k1 = so.shape(2) + 1;

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(2, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				so(i,j,k,seven_pt::ps) = 1.0*yh;
			}
		}
	}

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(2, i1)) {
				so(i,j,k,seven_pt::pw) = 1.0*xh;
			}
		}
	}

	for (auto k : range<len_t>(2, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				so(i,j,k,seven_pt::b) = 1.0*zh;
			}
		}
	}

	for (auto k : so.range(2)) {
		for (auto j : so.range(1)) {
			for (auto i : so.range(0)) {
				so(i,j,k,seven_pt::p) = 2.0*xh + 2.0*yh + 2.0*zh;
			}
		}
	}

	return so;
}

stencil_op<seven_pt> diag_diffusion(len_t nx, len_t ny, len_t nz,
                                    real_t dx, real_t dy, real_t dz)
{
	stencil_op<seven_pt> so(nx, ny, nz);

	so.set(0);

	real_t hx = 1.0/(so.len(0)-1);
	real_t hy = 1.0/(so.len(1)-1);
	real_t hz = 1.0/(so.len(2)-1);
	real_t xh=hy*hz/hx;
	real_t yh=hx*hz/hy;
	real_t zh=hx*hy/hz;
	len_t i1 = so.shape(0) + 1;
	len_t j1 = so.shape(1) + 1;
	len_t k1 = so.shape(2) + 1;

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(2, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				so(i,j,k,seven_pt::ps) = dy*yh;
			}
		}
	}

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(2, i1)) {
				so(i,j,k,seven_pt::pw) = dx*xh;
			}
		}
	}

	for (auto k : range<len_t>(2, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				so(i,j,k,seven_pt::b) = dz*zh;
			}
		}
	}

	for (auto k : so.range(2)) {
		for (auto j : so.range(1)) {
			for (auto i : so.range(0)) {
				so(i,j,k,seven_pt::p) = 2.0*dx*xh + 2.0*dy*yh + 2.0*dz*zh;
			}
		}
	}

	return so;
}


stencil_op<xxvii_pt> fe(len_t nx, len_t ny, len_t nz)
{
	stencil_op<xxvii_pt> so(nx, ny, nz);

	len_t i1 = so.shape(0) + 1;
	len_t j1 = so.shape(1) + 1;
	len_t k1 = so.shape(2) + 1;

	for (auto k : range<len_t>(1, k1)){
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(2, i1)) {
				so(i,j,k,xxvii_pt::pw) = 1.0;
			}
		}
	}

	for (auto k : range<len_t>(1, k1)){
		for (auto j : range<len_t>(2, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				so(i,j,k,xxvii_pt::ps) = 1.0;
			}
		}
	}

	for (auto k : range<len_t>(2, k1)){
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				so(i,j,k,xxvii_pt::b) = 1.0;
			}
		}
	}

	for (auto k : range<len_t>(1, k1)){
		for (auto j : range<len_t>(2, j1)) {
			for (auto i : range<len_t>(2, i1)) {
				so(i,j,k,xxvii_pt::pnw) = 1.0;
				so(i,j,k,xxvii_pt::psw) = 1.0;
			}
		}
	}

	for (auto k : range<len_t>(2, k1)){
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(2, i1)) {
				so(i,j,k,xxvii_pt::bw) = 1.0;
				so(i,j,k,xxvii_pt::be) = 1.0;
			}
		}
	}

	for (auto k : range<len_t>(2, k1)){
		for (auto j : range<len_t>(2, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				so(i,j,k,xxvii_pt::bn) = 1.0;
				so(i,j,k,xxvii_pt::bs) = 1.0;
			}
		}
	}

	for (auto k : range<len_t>(2, k1)){
		for (auto j : range<len_t>(2, j1)) {
			for (auto i : range<len_t>(2, i1)) {
				so(i,j,k,xxvii_pt::bnw) = 1.0;
				so(i,j,k,xxvii_pt::bne) = 1.0;
				so(i,j,k,xxvii_pt::bse) = 1.0;
				so(i,j,k,xxvii_pt::bsw) = 1.0;
			}
		}
	}

	for (auto k : so.range(2)) {
		for (auto j : so.range(1)) {
			for (auto i : so.range(0)) {
				so(i,j,k,xxvii_pt::p) = 26;
			}
		}
	}

	return so;
}

}}}
