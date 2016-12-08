#include <iostream>
#include <cmath>
#include <cstdlib>
#include <gtest/gtest.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include "pyrelax.h"

#include <boxmg/2d/gallery.h>
#include <boxmg/cycle/types.h>

TEST(relax_point, serial_2d) {
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	int nsweeps = 20;
	len_t nx = 12;
	len_t ny = nx;

	auto so = create_poisson(nx, ny);
	auto b = grid_func::zeros(nx, ny);
	auto x = grid_func::random(nx, ny);

	relax_stencil sor(nx, ny);

	auto kreg = so.get_registry();

	kreg->setup_relax(so, sor);

	for (auto i : range(nsweeps)) {
		(void)i;
		kreg->relax(so, x, b, sor, cycle::Dir::DOWN);
	}

	auto r = so.residual(x, b);
	auto rnorm = r.inf_norm();

	Py_Initialize();
	initpyrelax();

	double pynorm = gs_iter(nx, ny, nsweeps);

	// assert same order of magnitude
	ASSERT_EQ(static_cast<int>(std::log10(std::abs(rnorm))),
	          static_cast<int>(std::log10(std::abs(pynorm))));

	Py_Finalize();
}
