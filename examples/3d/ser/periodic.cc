#include <math.h>
#include <memory>
#include <array>
#include <iostream>

#include <cedar/types.h>
#include <cedar/3d/gallery.h>
#include <cedar/3d/solver.h>


static cedar::cdr3::stencil_op create_op(cedar::len_t nx, cedar::len_t ny, cedar::len_t nz,
                                         std::array<bool, 3> periodic)
{
	using namespace cedar;
	using namespace cedar::cdr3;

	auto so = stencil_op(nx, ny, nz);
	auto & sten = so.stencil();
	sten.five_pt() = true;

	sten.set(0);

	if (periodic[0]) nx--;
	if (periodic[1]) ny--;
	if (periodic[2]) nz--;
	real_t hx = 1.0/(nx+1);
	real_t hy = 1.0/(ny+1);
	real_t hz = 1.0/(nz+1);
	real_t xh=hy*hz/hx;
	real_t yh=hx*hz/hy;
	real_t zh=hx*hy/hz;
	len_t l = sten.shape(0);
	len_t m = sten.shape(1);
	len_t n = sten.shape(2);
	len_t i1 = sten.shape(0) + 1;
	len_t j1 = sten.shape(1) + 1;
	len_t k1 = sten.shape(2) + 1;
	len_t ibeg = 2;
	len_t jbeg = 2;
	len_t kbeg = 2;

	if (periodic[0]) ibeg--;
	if (periodic[1]) jbeg--;
	if (periodic[2]) kbeg--;

	auto & o = sten;

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(jbeg, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				o(i,j,k,dir::PS) = 1.0*yh;
			}
		}
	}

	for (auto k : range<len_t>(1, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(ibeg, i1)) {
				o(i,j,k,dir::PW) = 1.0*xh;
			}
		}
	}

	for (auto k : range<len_t>(kbeg, k1)) {
		for (auto j : range<len_t>(1, j1)) {
			for (auto i : range<len_t>(1, i1)) {
				o(i,j,k,dir::B) = 1.0*zh;
			}
		}
	}

	for (auto k : sten.range(2)) {
		for (auto j : sten.range(1)) {
			for (auto i : sten.range(0)) {
				o(i,j,k,dir::P) = 2.0*xh + 2.0*yh + 2.0*zh;
			}
		}
	}

	if (periodic[0]) {
		for (auto k : sten.range(2)) {
			for (auto j : sten.range(1)) {
				o(ibeg-1,j,k,dir::P ) = o(l,j,k,dir::P );
				o(ibeg-1,j,k,dir::PW) = o(l,j,k,dir::PW);
				o(ibeg-1,j,k,dir::PS) = o(l,j,k,dir::PS);
				o(ibeg-1,j,k,dir::B ) = o(l,j,k,dir::B );

				o(l+1,j,k,dir::P ) = o(ibeg,j,k,dir::P );
				o(l+1,j,k,dir::PW) = o(ibeg,j,k,dir::PW);
				o(l+1,j,k,dir::PS) = o(ibeg,j,k,dir::PS);
				o(l+1,j,k,dir::B ) = o(ibeg,j,k,dir::B );
			}
		}
	}

	if (periodic[1]) {
		for (auto k : sten.range(2)) {
			for (auto i : sten.range(0)) {
				o(i,jbeg-1,k,dir::P ) = o(i,m,k,dir::P );
				o(i,jbeg-1,k,dir::PW) = o(i,m,k,dir::PW);
				o(i,jbeg-1,k,dir::PS) = o(i,m,k,dir::PS);
				o(i,jbeg-1,k,dir::B ) = o(i,m,k,dir::B );

				o(i,m+1,k,dir::P ) = o(i,jbeg,k,dir::P );
				o(i,m+1,k,dir::PW) = o(i,jbeg,k,dir::PW);
				o(i,m+1,k,dir::PS) = o(i,jbeg,k,dir::PS);
				o(i,m+1,k,dir::B ) = o(i,jbeg,k,dir::B );
			}
		}
	}

	if (periodic[2]) {
		for (auto j: sten.range(1)) {
			for (auto i :sten.range(0)) {
				o(i,j,kbeg-1,dir::P ) = o(i,j,n,dir::P );
				o(i,j,kbeg-1,dir::PW) = o(i,j,n,dir::PW);
				o(i,j,kbeg-1,dir::PS) = o(i,j,n,dir::PS);
				o(i,j,kbeg-1,dir::B ) = o(i,j,n,dir::B );

				o(i,j,n+1,dir::P ) = o(i,j,kbeg,dir::P );
				o(i,j,n+1,dir::PW) = o(i,j,kbeg,dir::PW);
				o(i,j,n+1,dir::PS) = o(i,j,kbeg,dir::PS);
				o(i,j,n+1,dir::B ) = o(i,j,kbeg,dir::B );
			}
		}
	}

	return so;
}


static void set_problem(cedar::cdr3::grid_func & b, std::array<bool, 3> periodic)
{
	using namespace cedar;
	using namespace cedar::cdr3;

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


static void set_solution(cedar::cdr3::grid_func & q, std::array<bool, 3> periodic)
{
	using namespace cedar;

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


int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr3;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

	auto conf = std::make_shared<config::reader>();
	auto params = build_kernel_params(*conf);

	auto ndofs = conf->getvec<len_t>("grid.n");
	auto nx = ndofs[0];
	auto ny = ndofs[1];
	auto nz = ndofs[2];

	auto so = create_op(nx, ny, nz, params->periodic);
	grid_func b(nx, ny, nz);
	set_problem(b, params->periodic);

	solver bmg(std::move(so), conf);

	auto sol = bmg.solve(b);

	grid_func exact_sol(nx, ny, nz);

	set_solution(exact_sol, params->periodic);

	auto diff = exact_sol - sol;

	log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	log::status << "Finished test" << std::endl;

	MPI_Finalize();
	return 0;
}
