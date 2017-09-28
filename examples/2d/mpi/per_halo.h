#ifndef CEDAR_EXAMPLES_2D_PER_HALO_H
#define CEDAR_EXAMPLES_2D_PER_HALO_H

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/kernel/mpi/registry.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/tausch_exchanger.h>
#include <cedar/2d/mpi/msg_exchanger.h>
#include <cedar/kernel_params.h>

#include <cedar/util/time_log.h>

using namespace cedar;
using namespace cedar::cdr2;


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


static void draw_so(const cedar::cdr2::mpi::stencil_op<cedar::cdr2::five_pt> & so, std::string prefix)
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

	for (int i = 0; i < 3; ++i) {
		osv.emplace_back(std::make_shared<std::ofstream>("output/" + prefix + "-stencil-" +
		                                                 std::to_string(topo.coord(0)) +
		                                                 "." + std::to_string(topo.coord(1)) +
		                                                 "-" + std::to_string(i)));
	}

	for (auto j : so.grange(1)) {
		for (auto i : so.grange(0)) {
			*osv[0] << format(so(i,j,five_pt::c)) << " ";
			*osv[1] << format(so(i,j,five_pt::w)) << " ";
			*osv[2] << format(so(i,j,five_pt::s)) << " ";
		}
		for (int i = 0; i < 3; ++i) *(osv[i]) << '\n';
	}

	for (int i = 0; i < 3; ++i) osv[i]->close();
}


template<class halo_exchanger>
void run_test(config::reader & conf, std::shared_ptr<grid_topo> grid,
              mpi::stencil_op<five_pt> & so, mpi::grid_func & b)
{
	auto parms = build_kernel_params(conf);
	halo_exchanger halof(*parms, *grid);

	draw(b, "before");
	halof.exchange(b);
	MPI_Barrier(grid->comm);
	draw(b, "after");


	// draw_so(so, "before");
	// kreg.halo_stencil_exchange(so);
	// draw_so(so, "after");

	log::status << "Finished Test" << std::endl;

}

#endif
