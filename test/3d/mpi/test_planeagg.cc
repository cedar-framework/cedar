#include <gtest/gtest.h>

#include "line_agg.h"


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
