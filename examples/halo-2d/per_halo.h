#ifndef CEDAR_EXAMPLES_2D_PER_HALO_H
#define CEDAR_EXAMPLES_2D_PER_HALO_H

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/util/topo.h>
#include <cedar/kernel_params.h>

#include <cedar/util/time_log.h>

using namespace cedar;
using namespace cedar::cdr2;


static void fill_gfunc(cedar::cdr2::mpi::grid_func & b)
{
	auto & topo = b.grid();
	b.set(-1);
	for (auto j : b.range(1)) {
		for (auto i : b.range(0)) {
			b(i,j) = topo.coord(1)*topo.nproc(0) + topo.coord(0);
		}
	}
}


template<class sten>
static void fill_stencil(cedar::cdr2::mpi::stencil_op<sten> & so)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	so.set(-1);

	auto & topo = so.grid();

	for (auto j : so.range(1)) {
		for (auto i : so.range(0)) {
			for (int k = 0; k < stencil_ndirs<sten>::value; k++) {
				so(i,j,static_cast<sten>(k)) = 100*topo.coord(0) + 10*topo.coord(1) + k;
			}
		}
	}
}


static void draw(const cedar::cdr2::mpi::grid_func & b, std::string prefix)
{
	auto & topo = b.grid();
	std::ofstream os("output/" + prefix + "-gfunc-" + std::to_string(topo.coord(0)) +
	                 "." + std::to_string(topo.coord(1)));
	for (auto j : b.grange(1)) {
		for (auto i : b.grange(0)) {
			if (b(i,j) < 0)
				os << '*';
			else
				os << b(i,j);
			os << " ";
		}
		os << '\n';
	}

	os.close();
}


template<class sten>
static void draw_so(const cedar::cdr2::mpi::stencil_op<sten> & so, std::string prefix)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	std::vector<std::shared_ptr<std::ofstream>> osv;

	auto format = [](real_t v) -> std::string {
		if (v < 0.0)
			return "***";
		else {
			auto iv = static_cast<int>(v);
			std::string ret("");
			if (iv < 100)
				ret += "0";
			if (iv < 10)
				ret += "0";
			ret += std::to_string(static_cast<int>(v));
			return ret;
		}
	};

	auto & topo = so.grid();

	for (int i = 0; i < stencil_ndirs<sten>::value; ++i) {
		osv.emplace_back(std::make_shared<std::ofstream>("output/" + prefix + "-stencil-" +
		                                                 std::to_string(topo.coord(0)) +
		                                                 "." + std::to_string(topo.coord(1)) +
		                                                 "-" + std::to_string(i)));
	}

	for (auto j : so.grange(1)) {
		for (auto i : so.grange(0)) {
			for (int k = 0; k < stencil_ndirs<sten>::value; k++)
				*osv[0] << format(so(i,j,static_cast<sten>(k))) << " ";
		}
		for (int i = 0; i < stencil_ndirs<sten>::value; ++i) *(osv[i]) << '\n';
	}

	for (int i = 0; i < stencil_ndirs<sten>::value; ++i) osv[i]->close();
}

#endif
