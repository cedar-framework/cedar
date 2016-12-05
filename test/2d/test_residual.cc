#include <gtest/gtest.h>

#include <boxmg/2d/gallery.h>

TEST(residual, serial_2d) {
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	len_t nx = 5;
	len_t ny = nx;
	auto so = create_poisson(nx, ny);
	grid_func b(nx, ny);
	grid_func x(nx, ny);
	x.set(1);
	b.set(0);

	auto kreg = so.get_registry();

	kreg->set("residual", "fortran");
	auto r0 = so.residual(x, b);

	kreg->set("residual", "c++");
	auto r1 = so.residual(x, b);

	for (auto j : r0.grange(1)) {
		for (auto i : r0.grange(0)) {
			ASSERT_EQ(r0(i,j), r1(i,j));
		}
	}
}
