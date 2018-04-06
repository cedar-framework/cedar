#ifndef CEDAR_EXAMPLES_2D_PER_HALO_H
#define CEDAR_EXAMPLES_2D_PER_HALO_H

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/kernel/mpi/registry.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/tausch_exchanger.h>
#include <cedar/2d/mpi/msg_exchanger.h>
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


template<class halo_exchanger>
void run_test(config::reader & conf, std::shared_ptr<grid_topo> grid,
              mpi::stencil_op<five_pt> & so, mpi::grid_func & b)
{
	auto parms = build_kernel_params(conf);

	mpi::solver<five_pt, halo_exchanger> slv(so);
	auto kreg = slv.kernel_registry();

	draw(b, "before-0");
	draw_so(so, "before-0");
	kreg->halo_exchange(b);
	kreg->halo_stencil_exchange(so);
	draw(b, "after-0");
	draw_so(so, "after-0");


	for (std::size_t lvl = 1; lvl < slv.nlevels(); lvl++) {
		auto & level = slv.levels.get(lvl);
		fill_gfunc(level.b);
		fill_stencil(level.A);
		draw(level.b, "before-" + std::to_string(lvl));
		draw_so(level.A, "before-" + std::to_string(lvl));
		kreg->halo_exchange(level.b);
		kreg->halo_stencil_exchange(level.A);
		draw(level.b, "after-" + std::to_string(lvl));
		draw_so(level.A, "after-" + std::to_string(lvl));
	}

	log::status << "Finished Test" << std::endl;

}

#endif
