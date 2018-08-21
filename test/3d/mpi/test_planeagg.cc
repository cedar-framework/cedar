#include <gtest/gtest.h>

#include "line_agg.h"
#include "plane_agg.h"

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


template<>
cedar::cdr3::mpi::stencil_op<cedar::cdr3::seven_pt> create_op3(cedar::topo_ptr grid)
{
	return cedar::cdr3::mpi::gallery::poisson(grid);
}


template<>
cedar::cdr3::mpi::stencil_op<cedar::cdr3::xxvii_pt> create_op3(cedar::topo_ptr grid)
{
	return cedar::cdr3::mpi::gallery::fe(grid);
}


template<>
cedar::cdr2::mpi::stencil_op<cedar::cdr2::five_pt> create_op(cedar::topo_ptr grid)
{
	return cedar::cdr2::mpi::gallery::poisson(grid);
}


template<>
cedar::cdr2::mpi::stencil_op<cedar::cdr2::nine_pt> create_op(cedar::topo_ptr grid)
{
	return cedar::cdr2::mpi::gallery::fe(grid);
}


static void assert_tol(std::vector<cedar::real_t> & x0,
                       std::vector<cedar::real_t> & x1)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double tol = 1e-10;
	for (std::size_t i = 0; i < x0.size(); i++) {
		ASSERT_LT(std::abs(x0[i] - x1[i]), tol);
	}
}


TEST(MPIPlaneAgg, LineRelax5) {
	using namespace cedar;
	using namespace cedar::cdr2;

	auto xnoagg = line_agg<relax_dir::x, five_pt>(false, 31, 33, 11);
	auto xagg   = line_agg<relax_dir::x, five_pt>(true,  31, 33, 11);
	assert_tol(xnoagg, xagg);

	auto ynoagg = line_agg<relax_dir::y, five_pt>(false, 29, 28, 11);
	auto yagg   = line_agg<relax_dir::y, five_pt>(true,  29, 28, 11);
	assert_tol(ynoagg, yagg);
}


TEST(MPIPlaneAgg, LineRelax9) {
	using namespace cedar;
	using namespace cedar::cdr2;

	auto xnoagg = line_agg<relax_dir::x, nine_pt>(false, 31, 33, 10);
	auto xagg   = line_agg<relax_dir::x, nine_pt>(true,  31, 33, 10);
	assert_tol(xnoagg, xagg);

	auto ynoagg = line_agg<relax_dir::y, nine_pt>(false, 29, 28, 11);
	auto yagg   = line_agg<relax_dir::y, nine_pt>(true,  29, 28, 11);
	assert_tol(ynoagg, yagg);
}


TEST(MPIPlaneAgg, PlaneXY7) {
	using namespace cedar::cdr3;

	auto noagg = plane_agg<relax_dir::xy, seven_pt>(false, 24, 23, 13);
	auto agg   = plane_agg<relax_dir::xy, seven_pt>(true,  24, 23, 13);

	auto diff = noagg - agg;
	ASSERT_LT(std::abs(diff.inf_norm()), 1e-10);
}


TEST(MPIPlaneAgg, PlaneXZ7) {
	using namespace cedar::cdr3;

	auto noagg = plane_agg<relax_dir::xz, seven_pt>(false, 24, 21, 10);
	auto agg   = plane_agg<relax_dir::xz, seven_pt>(true,  24, 21, 10);

	auto diff = noagg - agg;
	ASSERT_LT(std::abs(diff.inf_norm()), 1e-10);
}


TEST(MPIPlaneAgg, PlaneYZ7) {
	using namespace cedar::cdr3;

	auto noagg = plane_agg<relax_dir::yz, seven_pt>(false, 14, 17, 11);
	auto agg   = plane_agg<relax_dir::yz, seven_pt>(true,  14, 17, 11);

	auto diff = noagg - agg;
	ASSERT_LT(std::abs(diff.inf_norm()), 1e-10);
}


TEST(MPIPlaneAgg, PlaneXY27) {
	using namespace cedar::cdr3;

	auto noagg = plane_agg<relax_dir::xy, xxvii_pt>(false, 25, 23, 13);
	auto agg   = plane_agg<relax_dir::xy, xxvii_pt>(true,  25, 23, 13);

	auto diff = noagg - agg;
	ASSERT_LT(std::abs(diff.inf_norm()), 1e-10);
}


TEST(MPIPlaneAgg, PlaneXZ27) {
	using namespace cedar::cdr3;

	auto noagg = plane_agg<relax_dir::xz, xxvii_pt>(false, 24, 21, 10);
	auto agg   = plane_agg<relax_dir::xz, xxvii_pt>(true,  24, 21, 10);

	auto diff = noagg - agg;
	ASSERT_LT(std::abs(diff.inf_norm()), 1e-10);
}


TEST(MPIPlaneAgg, PlaneYZ27) {
	using namespace cedar::cdr3;

	auto noagg = plane_agg<relax_dir::yz, xxvii_pt>(false, 14, 17, 11);
	auto agg   = plane_agg<relax_dir::yz, xxvii_pt>(true,  14, 17, 11);

	auto diff = noagg - agg;
	ASSERT_LT(std::abs(diff.inf_norm()), 1e-10);
}
