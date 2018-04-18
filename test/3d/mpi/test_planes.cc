#include "test_planes.h"

template<>
cedar::cdr3::stencil_op<cedar::cdr3::seven_pt> create_op_ser(cedar::len_t nx, cedar::len_t ny, cedar::len_t nz)
{
	return cedar::cdr3::gallery::poisson(nx, ny, nz);
}


template<>
cedar::cdr3::stencil_op<cedar::cdr3::xxvii_pt> create_op_ser(cedar::len_t nx, cedar::len_t ny, cedar::len_t nz)
{
	return cedar::cdr3::gallery::fe(nx, ny, nz);
}


template<>
cedar::cdr3::mpi::stencil_op<cedar::cdr3::seven_pt> create_op_mpi(cedar::topo_ptr grid)
{
	return cedar::cdr3::mpi::gallery::poisson(grid);
}


template<>
cedar::cdr3::mpi::stencil_op<cedar::cdr3::xxvii_pt> create_op_mpi(cedar::topo_ptr grid)
{
	return cedar::cdr3::mpi::gallery::fe(grid);
}


void set_rhs(cedar::cdr3::grid_func & b)
{
	for (auto k : b.range(2)) {
		for (auto j : b.range(1)) {
			for (auto i : b.range(0)) {
				b(i,j,k) = k*b.len(1)*b.len(0) + j * b.len(0) + i;
			}
		}
	}
}


void set_rhs(cedar::cdr3::mpi::grid_func & b)
{
	using namespace cedar;
	using namespace cedar::cdr3;
	auto & topo = b.grid();

	len_t igs = topo.is(0);
	len_t jgs = topo.is(1);
	len_t kgs = topo.is(2);

	len_t nx = topo.nglobal(0);
	len_t ny = topo.nglobal(1);

	for (auto k : b.range(2)) {
		for (auto j : b.range(1)) {
			for (auto i : b.range(0)) {
				len_t is = igs + i - 1;
				len_t js = jgs + j - 1;
				len_t ks = kgs + k - 1;
				b(i,j,k) = ks*nx*ny + js*nx + is;
			}
		}
	}
}


void set_sol(cedar::cdr3::grid_func & x)
{
	set_rhs(x);
	for (auto k : x.range(2)) {
		for (auto j : x.range(1)) {
			for (auto i : x.range(0)) {
				x(i,j,k) += 1;
			}
		}
	}
}


void set_sol(cedar::cdr3::mpi::grid_func & x)
{
	set_rhs(x);
	for (auto k : x.range(2)) {
		for (auto j : x.range(1)) {
			for (auto i : x.range(0)) {
				x(i,j,k) += 1;
			}
		}
	}
}


TEST(MPIPlanes, XY7) {
	using namespace cedar;
	using namespace cedar::cdr3;

	plane_test<relax_dir::xy, seven_pt>(50, 51, 5);
}


TEST(MPIPlanes, XY27) {
	using namespace cedar;
	using namespace cedar::cdr3;

	plane_test<relax_dir::xy, xxvii_pt>(50, 51, 6);
}


TEST(MPIPlanes, XZ7) {
	using namespace cedar;
	using namespace cedar::cdr3;

	plane_test<relax_dir::xz, seven_pt>(50, 5, 51);
}


TEST(MPIPlanes, XZ27) {
	using namespace cedar;
	using namespace cedar::cdr3;

	plane_test<relax_dir::xz, xxvii_pt>(50, 6, 49);
}


TEST(MPIPlanes, YZ7) {
	using namespace cedar;
	using namespace cedar::cdr3;

	plane_test<relax_dir::yz, seven_pt>(4, 49, 51);
}


TEST(MPIPlanes, YZ27) {
	using namespace cedar;
	using namespace cedar::cdr3;

	plane_test<relax_dir::yz, xxvii_pt>(5, 50, 48);
}
