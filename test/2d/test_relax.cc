#include <iostream>
#include <cmath>
#include <cstdlib>
#include <gtest/gtest.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include "pyrelax.h"

#include <boxmg/2d/gallery.h>
#include <boxmg/cycle/types.h>


TEST(relax_point, serial_2d_fivept) {
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	int nsweeps = 7;
	len_t nx = 31;
	len_t ny = nx;

	auto so = create_poisson(nx, ny);
	auto b = grid_func::zeros(nx, ny);
	auto x = grid_func::ones(nx, ny);

	relax_stencil sor(nx, ny);

	auto kreg = so.get_registry();

	kreg->setup_relax(so, sor);

	for (auto i : range(nsweeps)) {
		(void)i;
		kreg->relax(so, x, b, sor, cycle::Dir::DOWN);
	}

	for (auto i : range(nsweeps)) {
		(void)i;
		kreg->relax(so, x, b, sor, cycle::Dir::UP);
	}

	grid_func pyx(nx, ny, 0);

	gs_iter(nx, ny, nsweeps, 5, pyx.data());

	auto diff = x - pyx;

	auto tol = 1e-10;

	ASSERT_LT(std::abs(diff.inf_norm()), tol);
}


TEST(relax_point, serial_2d_ninept) {
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	int nsweeps = 3;
	len_t nx = 37;
	len_t ny = nx;

	auto so = create_fe(nx, ny);
	auto b = grid_func::zeros(nx, ny);
	auto x = grid_func::ones(nx, ny);

	relax_stencil sor(nx, ny);

	auto kreg = so.get_registry();

	kreg->setup_relax(so, sor);

	for (auto i : range(nsweeps)) {
		(void)i;
		kreg->relax(so, x, b, sor, cycle::Dir::DOWN);
	}


	for (auto i : range(nsweeps)) {
		(void)i;
		kreg->relax(so, x, b, sor, cycle::Dir::UP);
	}

	grid_func pyx(nx, ny, 0);
	gs_iter(nx, ny, nsweeps, 9, pyx.data());


	auto diff = x - pyx;

	auto tol = 1e-10;

	ASSERT_LT(std::abs(diff.inf_norm()), tol);
}


TEST(relax_line_x, serial_2d_fivept) {
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	int nsweeps = 15;
	len_t nx = 132;
	len_t ny = 132;

	auto so = create_diag_diffusion(nx, ny, 1, .0001);
	auto b = grid_func::zeros(nx, ny);
	auto x = grid_func::random(nx, ny);

	auto & sten = so.stencil();

	relax_stencil sor(nx, ny);

	auto kreg = so.get_registry();

	kreg->setup_relax_x(so, sor);

	grid_func tmp(nx, ny);

	for (auto i : range(nsweeps)) {
		(void)i;
		kreg->relax_lines_x(so, x, b, sor, tmp, cycle::Dir::DOWN);
	}

	double max_high_freq = get_high_freq(nx, ny, x.data());

	auto tol = 1e-10;
	ASSERT_LT(max_high_freq, tol);
}


TEST(relax_line_y, serial_2d_fivept) {
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	int nsweeps = 15;
	len_t nx = 132;
	len_t ny = 132;

	auto so = create_diag_diffusion(nx, ny, .0001, 1);
	auto b = grid_func::zeros(nx, ny);
	auto x = grid_func::random(nx, ny);

	auto & sten = so.stencil();

	relax_stencil sor(nx, ny);

	auto kreg = so.get_registry();

	kreg->setup_relax_y(so, sor);

	grid_func tmp(nx, ny);

	for (auto i : range(nsweeps)) {
		(void)i;
		kreg->relax_lines_y(so, x, b, sor, tmp, cycle::Dir::DOWN);
	}

	double max_high_freq = get_high_freq(nx, ny, x.data());

	auto tol = 1e-10;
	ASSERT_LT(max_high_freq, tol);
}


int main(int argc, char *argv[])
{
	Py_Initialize();
	initpyrelax();

	::testing::InitGoogleTest(&argc, argv);
	auto ret = RUN_ALL_TESTS();

	Py_Finalize();
	return ret;
}
