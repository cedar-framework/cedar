#include <iostream>
#include <cmath>
#include <cstdlib>
#include <gtest/gtest.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include "pyrelax.h"

#include <cedar/2d/gallery.h>
#include <cedar/cycle/types.h>
#include <cedar/2d/kernel_manager.h>


TEST(SerialRelax2, Point5) {
	using namespace cedar;
	using namespace cedar::cdr2;

	int nsweeps = 7;
	len_t nx = 31;
	len_t ny = nx;

	config conf("");
	log::init(conf);

	auto kman = build_kernel_manager(conf);

	auto so = gallery::poisson(nx, ny);
	auto b = grid_func::zeros(nx, ny);
	auto x = grid_func::ones(nx, ny);

	relax_stencil sor(nx, ny);

	kman->setup<point_relax>(so, sor);

	for (auto i : range(nsweeps)) {
		(void)i;
		kman->run<point_relax>(so, x, b, sor, cycle::Dir::DOWN);
	}

	for (auto i : range(nsweeps)) {
		(void)i;
		kman->run<point_relax>(so, x, b, sor, cycle::Dir::UP);
	}

	grid_func pyx(nx, ny, 0);

	gs_iter(nx, ny, nsweeps, 5, pyx.data());

	auto diff = x - pyx;

	auto tol = 1e-10;

	ASSERT_LT(std::abs(diff.inf_norm()), tol);
}


TEST(SerialRelax2, Point9) {
	using namespace cedar;
	using namespace cedar::cdr2;

	int nsweeps = 3;
	len_t nx = 37;
	len_t ny = nx;

	auto so = gallery::fe(nx, ny);
	auto b = grid_func::zeros(nx, ny);
	auto x = grid_func::ones(nx, ny);

	relax_stencil sor(nx, ny);

	config conf("");
	log::init(conf);
	auto kman = build_kernel_manager(conf);

	kman->setup<point_relax>(so, sor);

	for (auto i : range(nsweeps)) {
		(void)i;
		kman->run<point_relax>(so, x, b, sor, cycle::Dir::DOWN);
	}


	for (auto i : range(nsweeps)) {
		(void)i;
		kman->run<point_relax>(so, x, b, sor, cycle::Dir::UP);
	}

	grid_func pyx(nx, ny, 0);
	gs_iter(nx, ny, nsweeps, 9, pyx.data());


	auto diff = x - pyx;

	auto tol = 1e-10;

	ASSERT_LT(std::abs(diff.inf_norm()), tol);
}


TEST(SerialRelax2, LineX5) {
	using namespace cedar;
	using namespace cedar::cdr2;

	int nsweeps = 15;
	len_t nx = 132;
	len_t ny = 132;

	config conf("");
	log::init(conf);
	auto kman = build_kernel_manager(conf);

	auto so = gallery::diag_diffusion(nx, ny, 1, .0001);
	auto b = grid_func::zeros(nx, ny);
	auto x = grid_func::random(nx, ny);

	relax_stencil sor(nx, ny);

	kman->setup<line_relax<relax_dir::x>>(so, sor);

	grid_func tmp(nx, ny);

	for (auto i : range(nsweeps)) {
		(void)i;
		kman->run<line_relax<relax_dir::x>>(so, x, b, sor, tmp, cycle::Dir::DOWN);
	}

	double max_high_freq = get_high_freq(nx, ny, x.data());

	auto tol = 1e-10;
	ASSERT_LT(max_high_freq, tol);
}


TEST(SerialRelax2, LineY5) {
	using namespace cedar;
	using namespace cedar::cdr2;

	int nsweeps = 15;
	len_t nx = 132;
	len_t ny = 132;

	config conf("");
	log::init(conf);
	auto kman = build_kernel_manager(conf);

	auto so = gallery::diag_diffusion(nx, ny, .0001, 1);
	auto b = grid_func::zeros(nx, ny);
	auto x = grid_func::random(nx, ny);

	relax_stencil sor(nx, ny);

	kman->setup<line_relax<relax_dir::y>>(so, sor);

	grid_func tmp(nx, ny);

	for (auto i : range(nsweeps)) {
		(void)i;
		kman->run<line_relax<relax_dir::y>>(so, x, b, sor, tmp, cycle::Dir::DOWN);
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
