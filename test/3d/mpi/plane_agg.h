#ifndef CEDAR_TEST_PLANEAGG_H
#define CEDAR_TEST_PLANEAGG_H

#include <cedar/3d/mpi/gallery.h>
#include <cedar/3d/util/topo.h>
#include <cedar/3d/mpi/types.h>
#include <cedar/3d/mpi/kernel_manager.h>


using namespace cedar;
using namespace cedar::cdr3;


void set_rhs(cdr3::mpi::grid_func & b);
void set_sol(cdr3::mpi::grid_func & x);

template<class sten>
cdr3::mpi::stencil_op<sten> create_op3(topo_ptr grid);


namespace cedar { namespace cdr3 {

template<relax_dir rdir, class sten>
mpi::grid_func plane_agg(bool aggregate, len_t nx, len_t ny, len_t nz)
{
	std::size_t nsweeps = 3;

	config conf("test-planeagg.json");
	log::init(conf);
	conf.set("solver.plane-agg", aggregate);

	auto kman = mpi::build_kernel_manager(conf);

	auto grid = util::create_topo(MPI_COMM_WORLD, nx, ny, nz);
	auto so = create_op3<sten>(grid);
	mpi::grid_func x(grid), b(grid);
	set_rhs(b);
	set_sol(x);

	{ // setup halo
		auto & halo_service = kman->services().get<mpi::halo_exchange>();
		halo_service.setup({{grid}});
	}

	kman->setup<mpi::plane_relax<rdir>>(so);

	for (auto i : range<std::size_t>(nsweeps)) {
		kman->run<mpi::plane_relax<rdir>>(so, x, b, cycle::Dir::DOWN);
		kman->run<mpi::plane_relax<rdir>>(so, x, b, cycle::Dir::UP);
	}

	return x;
}

}}
#endif
