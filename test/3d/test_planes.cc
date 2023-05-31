#include <gtest/gtest.h>
#include <cedar/3d/gallery.h>
#include <cedar/cycle/types.h>
#include <cedar/3d/types.h>
#include <cedar/3d/kernel_manager.h>

#include <cedar/3d/relax_planes.h>

#include <Python.h>
#include <numpy/arrayobject.h>
#include "pyplanes.h"

TEST(SerialPlanes, XY7) {
	using namespace cedar;
	using namespace cedar::cdr3;

	len_t nx = 31;
	len_t ny = 35;
	len_t nz = 5;

	config conf("test-planes-ser.json");
	log::init(conf);

	auto kman = build_kernel_manager(conf);

	auto so = gallery::poisson(nx, ny, nz);
	grid_func x(nx,ny,nz), b(nx,ny,nz), pyx(nx,ny,nz);
	cdr2::grid_func pyb(nx, ny);

	for (auto k : x.range(2)) {
		for (auto j : x.range(1)) {
			for (auto i : x.range(0)) {
				x(i,j,k) = k*nx*ny + j*nx + i;
				b(i,j,k) = k*nx*ny + j*nx + i;
			}
		}
	}

	kman->setup<plane_relax<relax_dir::xy>>(so);
	kman->run<plane_relax<relax_dir::xy>>(so, x, b, cycle::Dir::DOWN);

	float tol = 1e-8;
	for (auto k : range<len_t>(2, nz + 1, 2)) {
		copy_rhs<relax_dir::xy>(so, x, b, pyb, k);
		pyplane(nx, ny, nz, pyx.data(), pyb.data(), k);
		for (auto j : x.range(1)) {
			for (auto i : x.range(0)) {
				ASSERT_LT(std::abs(x(i,j,k) - pyx(i,j,k)), tol);
			}
		}
	}
}


TEST(SerialPlanes, XY27) {
	using namespace cedar;
	using namespace cedar::cdr3;

	len_t nx = 31;
	len_t ny = 35;
	len_t nz = 5;

	config conf("test-planes-ser.json");
	log::init(conf);

	auto kman = build_kernel_manager(conf);

	auto so = gallery::fe(nx, ny, nz);
	grid_func x(nx,ny,nz), b(nx,ny,nz), pyx(nx,ny,nz);
	cdr2::grid_func pyb(nx, ny);

	for (auto k : x.range(2)) {
		for (auto j : x.range(1)) {
			for (auto i : x.range(0)) {
				x(i,j,k) = k*nx*ny + j*nx + i;
				b(i,j,k) = k*nx*ny + j*nx + i;
			}
		}
	}

	kman->setup<plane_relax<relax_dir::xy>>(so);
	kman->run<plane_relax<relax_dir::xy>>(so, x, b, cycle::Dir::DOWN);

	float tol = 1e-8;
	for (auto k : range<len_t>(2, nz + 1, 2)) {
		copy_rhs<relax_dir::xy>(so, x, b, pyb, k);
		pyplane27(nx, ny, nz, pyx.data(), pyb.data(), k);
		for (auto j : x.range(1)) {
			for (auto i : x.range(0)) {
				ASSERT_LT(std::abs(x(i,j,k) - pyx(i,j,k)), tol);
			}
		}
	}
}


TEST(SerialPlanes, XZ7) {
	using namespace cedar;
	using namespace cedar::cdr3;

	len_t nx = 31;
	len_t ny = 5;
	len_t nz = 35;

	config conf("test-planes-ser.json");
	log::init(conf);

	auto kman = build_kernel_manager(conf);

	auto so = gallery::poisson(nx, ny, nz);
	grid_func x(nx,ny,nz), b(nx,ny,nz), pyx(nx,ny,nz);
	cdr2::grid_func pyb(nx, nz);

	for (auto k : x.range(2)) {
		for (auto j : x.range(1)) {
			for (auto i : x.range(0)) {
				x(i,j,k) = k*nx*ny + j*nx + i;
				b(i,j,k) = k*nx*ny + j*nx + i;
			}
		}
	}

	kman->setup<plane_relax<relax_dir::xz>>(so);
	kman->run<plane_relax<relax_dir::xz>>(so, x, b, cycle::Dir::DOWN);

	float tol = 1e-8;
	for (auto j : range<len_t>(2, ny + 1, 2)) {
		copy_rhs<relax_dir::xz>(so, x, b, pyb, j);
		pyplane_xz(nx, ny, nz, pyx.data(), pyb.data(), j);
		for (auto k : x.range(2)) {
			for (auto i : x.range(0)) {
				ASSERT_LT(std::abs(x(i,j,k) - pyx(i,j,k)), tol);
			}
		}
	}
}


TEST(SerialPlanes, XZ27) {
	using namespace cedar;
	using namespace cedar::cdr3;

	len_t nx = 31;
	len_t ny = 5;
	len_t nz = 35;

	config conf("test-planes-ser.json");
	log::init(conf);

	auto kman = build_kernel_manager(conf);

	auto so = gallery::fe(nx, ny, nz);
	grid_func x(nx,ny,nz), b(nx,ny,nz), pyx(nx,ny,nz);
	cdr2::grid_func pyb(nx, nz);

	for (auto k : x.range(2)) {
		for (auto j : x.range(1)) {
			for (auto i : x.range(0)) {
				x(i,j,k) = k*nx*ny + j*nx + i;
				b(i,j,k) = k*nx*ny + j*nx + i;
			}
		}
	}

	kman->setup<plane_relax<relax_dir::xz>>(so);
	kman->run<plane_relax<relax_dir::xz>>(so, x, b, cycle::Dir::DOWN);

	float tol = 1e-8;
	for (auto j : range<len_t>(2, ny + 1, 2)) {
		copy_rhs<relax_dir::xz>(so, x, b, pyb, j);
		pyplane_xz27(nx, ny, nz, pyx.data(), pyb.data(), j);
		for (auto k : x.range(2)) {
			for (auto i : x.range(0)) {
				ASSERT_LT(std::abs(x(i,j,k) - pyx(i,j,k)), tol);
			}
		}
	}
}


TEST(SerialPlanes, YZ7) {
	using namespace cedar;
	using namespace cedar::cdr3;

	len_t nx = 5;
	len_t ny = 31;
	len_t nz = 35;

	config conf("test-planes-ser.json");
	log::init(conf);

	auto kman = build_kernel_manager(conf);

	auto so = gallery::poisson(nx, ny, nz);
	grid_func x(nx,ny,nz), b(nx,ny,nz), pyx(nx,ny,nz);
	cdr2::grid_func pyb(ny, nz);

	for (auto k : x.range(2)) {
		for (auto j : x.range(1)) {
			for (auto i : x.range(0)) {
				x(i,j,k) = k*nx*ny + j*nx + i;
				b(i,j,k) = k*nx*ny + j*nx + i;
			}
		}
	}

	kman->setup<plane_relax<relax_dir::yz>>(so);
	kman->run<plane_relax<relax_dir::yz>>(so, x, b, cycle::Dir::DOWN);

	float tol = 1e-8;
	for (auto i : range<len_t>(2, nx + 1, 2)) {
		copy_rhs<relax_dir::yz>(so, x, b, pyb, i);
		pyplane_yz(nx, ny, nz, pyx.data(), pyb.data(), i);
		for (auto k : x.range(2)) {
			for (auto j : x.range(1)) {
				ASSERT_LT(std::abs(x(i,j,k) - pyx(i,j,k)), tol);
			}
		}
	}
}


TEST(SerialPlanes, YZ27) {
	using namespace cedar;
	using namespace cedar::cdr3;

	len_t nx = 5;
	len_t ny = 31;
	len_t nz = 35;

	config conf("test-planes-ser.json");
	log::init(conf);

	auto kman = build_kernel_manager(conf);

	auto so = gallery::fe(nx, ny, nz);
	grid_func x(nx,ny,nz), b(nx,ny,nz), pyx(nx,ny,nz);
	cdr2::grid_func pyb(ny, nz);

	for (auto k : x.range(2)) {
		for (auto j : x.range(1)) {
			for (auto i : x.range(0)) {
				x(i,j,k) = k*nx*ny + j*nx + i;
				b(i,j,k) = k*nx*ny + j*nx + i;
			}
		}
	}

	kman->setup<plane_relax<relax_dir::yz>>(so);
	kman->run<plane_relax<relax_dir::yz>>(so, x, b, cycle::Dir::DOWN);

	float tol = 1e-8;
	for (auto i : range<len_t>(2, nx + 1, 2)) {
		copy_rhs<relax_dir::yz>(so, x, b, pyb, i);
		pyplane_yz27(nx, ny, nz, pyx.data(), pyb.data(), i);
		for (auto k : x.range(2)) {
			for (auto j : x.range(1)) {
				ASSERT_LT(std::abs(x(i,j,k) - pyx(i,j,k)), tol);
			}
		}
	}
}


int main(int argc, char *argv[])
{
	PyImport_AppendInittab("pyplanes", PyInit_pyplanes);
	Py_Initialize();
	PyInit_pyplanes();
	PyImport_ImportModule("pyplanes");

	::testing::InitGoogleTest(&argc, argv);
	auto ret = RUN_ALL_TESTS();

	Py_Finalize();
	return ret;
}
