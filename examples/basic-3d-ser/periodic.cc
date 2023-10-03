#include <math.h>
#include <memory>
#include <array>
#include <iostream>

#include <cedar/types.h>
#include <cedar/3d/gallery.h>
#include <cedar/3d/solver.h>

using namespace cedar;
using namespace cedar::cdr3;

static stencil_op<seven_pt> create_op(len_t nx, len_t ny, len_t nz,
                                      std::array<bool, 3> periodic)
{
	stencil_op<seven_pt> so(nx, ny, nz);

	so.set(0);

	if (periodic[0]) nx--;
	if (periodic[1]) ny--;
	if (periodic[2]) nz--;
	real_t hx = 1.0/(nx+1);
	real_t hy = 1.0/(ny+1);
	real_t hz = 1.0/(nz+1);
	real_t xh=hy*hz/hx;
	real_t yh=hx*hz/hy;
	real_t zh=hx*hy/hz;
	len_t l = so.shape(0);
	len_t m = so.shape(1);
	len_t n = so.shape(2);
	len_t i1 = so.shape(0) + 1;
	len_t j1 = so.shape(1) + 1;
	len_t k1 = so.shape(2) + 1;
	len_t ibeg = 2;
	len_t jbeg = 2;
	len_t kbeg = 2;

	if (periodic[0]) ibeg--;
	if (periodic[1]) jbeg--;
	if (periodic[2]) kbeg--;

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(jbeg, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				so(i,j,k,seven_pt::ps) = 1.0*yh;
			}
		}
	}

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(ibeg, i1)) {
				so(i,j,k,seven_pt::pw) = 1.0*xh;
			}
		}
	}

	for (auto k : range<len_t>(kbeg, k1)) {
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

	if (periodic[0]) {
		for (auto k : so.grange(2)) {
			for (auto j : so.grange(1)) {
				so(ibeg-1,j,k,seven_pt::p ) = so(l,j,k,seven_pt::p );
				so(ibeg-1,j,k,seven_pt::pw) = so(l,j,k,seven_pt::pw);
				so(ibeg-1,j,k,seven_pt::ps) = so(l,j,k,seven_pt::ps);
				so(ibeg-1,j,k,seven_pt::b ) = so(l,j,k,seven_pt::b );

				so(l+1,j,k,seven_pt::p ) = so(ibeg,j,k,seven_pt::p );
				so(l+1,j,k,seven_pt::pw) = so(ibeg,j,k,seven_pt::pw);
				so(l+1,j,k,seven_pt::ps) = so(ibeg,j,k,seven_pt::ps);
				so(l+1,j,k,seven_pt::b ) = so(ibeg,j,k,seven_pt::b );
			}
		}
	}

	if (periodic[1]) {
		for (auto k : so.grange(2)) {
			for (auto i : so.grange(0)) {
				so(i,jbeg-1,k,seven_pt::p ) = so(i,m,k,seven_pt::p );
				so(i,jbeg-1,k,seven_pt::pw) = so(i,m,k,seven_pt::pw);
				so(i,jbeg-1,k,seven_pt::ps) = so(i,m,k,seven_pt::ps);
				so(i,jbeg-1,k,seven_pt::b ) = so(i,m,k,seven_pt::b );

				so(i,m+1,k,seven_pt::p ) = so(i,jbeg,k,seven_pt::p );
				so(i,m+1,k,seven_pt::pw) = so(i,jbeg,k,seven_pt::pw);
				so(i,m+1,k,seven_pt::ps) = so(i,jbeg,k,seven_pt::ps);
				so(i,m+1,k,seven_pt::b ) = so(i,jbeg,k,seven_pt::b );
			}
		}
	}

	if (periodic[2]) {
		for (auto j: so.grange(1)) {
			for (auto i :so.grange(0)) {
				so(i,j,kbeg-1,seven_pt::p ) = so(i,j,n,seven_pt::p );
				so(i,j,kbeg-1,seven_pt::pw) = so(i,j,n,seven_pt::pw);
				so(i,j,kbeg-1,seven_pt::ps) = so(i,j,n,seven_pt::ps);
				so(i,j,kbeg-1,seven_pt::b ) = so(i,j,n,seven_pt::b );

				so(i,j,n+1,seven_pt::p ) = so(i,j,kbeg,seven_pt::p );
				so(i,j,n+1,seven_pt::pw) = so(i,j,kbeg,seven_pt::pw);
				so(i,j,n+1,seven_pt::ps) = so(i,j,kbeg,seven_pt::ps);
				so(i,j,n+1,seven_pt::b ) = so(i,j,kbeg,seven_pt::b );
			}
		}
	}

	return so;
}


static void set_problem(grid_func & b, std::array<bool, 3> periodic)
{
	const double pi = M_PI;

	auto rhs = [pi](real_t x, real_t y, real_t z) {
		return 12*(pi*pi)*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z);
	};

	b.set(0);

	len_t nx = b.len(0) - 2;
	len_t ny = b.len(1) - 2;
	len_t nz = b.len(2) - 2;

	if (periodic[0]) nx--;
	if (periodic[1]) ny--;
	if (periodic[2]) nz--;

	real_t hx = 1.0 / (nx + 1);
	real_t hy = 1.0 / (ny + 1);
	real_t hz = 1.0 / (nz + 1);

	real_t h2 = hx*hy*hz;

	for (auto k : b.range(2)) {
		for (auto j : b.range(1)) {
			for (auto i : b.range(0)) {
				real_t x = i*hx;
				real_t y = j*hy;
				real_t z = k*hz;

				b(i,j,k) = rhs(x,y,z) * h2;
			}
		}
	}

	if (periodic[0]) {
		for (auto k : b.grange(2)) {
			for (auto j : b.grange(1)) {
				b(0           ,j,k) = b(b.shape(0),j,k);
				b(b.shape(0)+1,j,k) = b(1         ,j,k);
			}
		}
	}

	if (periodic[1]) {
		for (auto k : b.grange(2)) {
			for (auto i : b.grange(0)) {
				b(i,0,           k) = b(i,b.shape(1),k);
				b(i,b.shape(1)+1,k) = b(i,1         ,k);
			}
		}
	}

	if (periodic[2]) {
		for (auto j : b.grange(1)) {
			for (auto i : b.grange(0)) {
				b(i,j,0           ) = b(i,j,b.shape(2));
				b(i,j,b.shape(2)+1) = b(i,j,1         );
			}
		}
	}
}


static void set_solution(grid_func & q, std::array<bool, 3> periodic)
{
	const double pi = M_PI;

	auto sol = [pi](real_t x, real_t y, real_t z) {
		return sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z);
	};

	len_t nx = q.len(0) - 2;
	len_t ny = q.len(1) - 2;
	len_t nz = q.len(2) - 2;
	if (periodic[0]) nx--;
	if (periodic[1]) ny--;
	if (periodic[2]) nz--;

	real_t hx = 1.0 / (nx + 1);
	real_t hy = 1.0 / (ny + 1);
	real_t hz = 1.0 / (nz + 1);

	for (auto k : q.range(2)) {
		for (auto j : q.range(1)) {
			for (auto i : q.range(0)) {
				real_t x = i*hx;
				real_t y = j*hy;
				real_t z = k*hz;

				q(i,j,k) = sol(x,y,z);
			}
		}
	}
}


int main()
{
	auto conf = std::make_shared<config>();
	auto params = build_kernel_params(*conf);

	auto ndofs = conf->getvec<len_t>("grid.n");
	auto nx = ndofs[0];
	auto ny = ndofs[1];
	auto nz = ndofs[2];

	auto so = create_op(nx, ny, nz, params->periodic);
	grid_func b(nx, ny, nz);
	set_problem(b, params->periodic);

	solver<seven_pt> bmg(so, conf);

	{
		std::ofstream ffile("output/ser-fine");
		std::ofstream rfile("output/ser-restrict");
		std::ofstream cfile("output/ser-coarse");

		ffile << bmg.levels.get<seven_pt>(0).A;
		rfile << bmg.levels.get(1).P;
		cfile << bmg.levels.get(1).A;

		ffile.close();
		rfile.close();
		cfile.close();
	}

	auto sol = bmg.solve(b);

	grid_func exact_sol(nx, ny, nz);

	set_solution(exact_sol, params->periodic);

	auto diff = exact_sol - sol;

	log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	log::status << "Finished test" << std::endl;

	return 0;
}
