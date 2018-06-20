#include <mpi.h>
#include <omp.h>

#include <iostream>

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/mpi/gallery.h>
#include <cedar/2d/mpi/kernel_manager.h>

#include <cedar/util/time_log.h>

#include "exchangers.h"


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


static void fill_solution(std::vector<cedar::cdr2::mpi::grid_func> & xs)
{
	auto nplanes = xs.size();
	auto grid = xs[0].grid_ptr();

	for (auto k : cedar::range(nplanes)) {
		for (auto j : xs[k].range(1)) {
			for (auto i : xs[k].range(0)) {
				xs[k](i,j) = k * grid->nproc(0)*grid->nproc(1) + grid->coord(1) * grid->nproc(0) + grid->coord(0) + 1;
			}
		}
	}
}


int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr2;

	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);

	timer_init(MPI_COMM_WORLD);

	auto conf = std::make_shared<config>("config.json");
	log::init(*conf);

	auto nplanes = conf->get<int>("nplanes");
	auto nsweeps = conf->get<int>("nsweeps");
	auto aggregate = conf->get<bool>("aggregate");

	auto grid = util::create_topo(*conf);
	auto so = mpi::gallery::fe(grid);
	relax_stencil sor(grid->nlocal(0) - 2, grid->nlocal(1) - 2);
	mpi::grid_func b(grid);
	b.set(1.0);
	std::vector<mpi::grid_func> xs;
	std::vector<mpi::kman_ptr> kmans;
	kmans.reserve(nplanes);
	xs.reserve(nplanes);
	real_t *x_data = new real_t[grid->nlocal(0)*grid->nlocal(1)*nplanes];
	real_t *x_curr = x_data;
	for (auto i : range(nplanes)) {
		(void)i;
		xs.emplace_back(x_curr, grid);
		xs.back().set(-1);
		x_curr += grid->nlocal(0) * grid->nlocal(1);

		kmans.push_back(mpi::build_kernel_manager(*conf));

		auto kman = kmans.back();
		kman->add<mpi::halo_exchange, mpi::noop_exchange>("noop");
		kman->add<mpi::halo_exchange, mpi::master_exchange>("master", nplanes);
	}

	if (aggregate) {
		for (auto i : range(1, nplanes)) {
			kmans[i]->set<mpi::halo_exchange>("noop");
		}
		kmans[0]->set<mpi::halo_exchange>("master");

		for (auto kman : kmans) {
			auto hex = kman->get_ptr<mpi::halo_exchange>();
			kman->add_halo(hex.get());
		}
	}

	{ // setup halo
		std::vector<topo_ptr> topos{{grid}};
		kmans[0]->setup<mpi::halo_exchange>(topos);
		if (not aggregate) {
			for (auto i : range(1, nplanes))
				kmans[i]->setup<mpi::halo_exchange>(topos);
		}
	}

	if (aggregate) {
		kmans.back()->setup<mpi::point_relax>(so, sor);

		timer_begin("relax");
		#pragma omp parallel num_threads(nplanes)
		{
			int tid = omp_get_thread_num();
			for (auto sweep : range(nsweeps)) {
				kmans[tid]->run<mpi::point_relax>(so, xs[tid], b, sor, cycle::Dir::DOWN);
			}
			// kmans[tid]->run<mpi::halo_exchange>(xs[tid]);
		}
		timer_end("relax");
	} else {
		for (auto i : range(nplanes)) {
			kmans[i]->setup<mpi::point_relax>(so, sor);
		}

		timer_begin("relax");
		for (auto sweep : range(nsweeps)) {
			for (auto i : range(nplanes)) {
				kmans[i]->run<mpi::point_relax>(so, xs[i], b, sor, cycle::Dir::DOWN);
			}
		}
		timer_end("relax");
	}

	timer_save("thread-times.json");

	delete[] x_data;

	MPI_Finalize();
	return 0;
}
