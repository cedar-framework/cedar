#ifndef CEDAR_3D_MPI_RELAX_PLANES_H
#define CEDAR_3D_MPI_RELAX_PLANES_H

#include <tuple>
#include <type_traits>

#include <cedar/types.h>
#include <cedar/kernels/plane_relax.h>
#include <cedar/3d/mpi/types.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/3d/mpi/kernel_manager.h>
#include <cedar/3d/mpi/plane_ult.h>
#include <cedar/3d/mpi/plane_team.h>
#include <cedar/3d/mpi/plane_util.h>

namespace cedar { namespace cdr3 { namespace mpi {

template<relax_dir rdir>
class planes : public kernels::plane_relax<stypes, rdir>
{
public:
	planes() : aggregate(false) {}

	void setup(const stencil_op<seven_pt> & so) override
	{
		level_teams.emplace_back();
		timer_pause();
		this->setup_impl(so, fine_planes, fine_threads, level_teams.back());
		timer_play();
	}
	void setup(const stencil_op<xxvii_pt> & so) override
	{
		const auto & grid = so.grid();
		level_threads.reserve(grid.nlevel()-1);
		level_teams.reserve(grid.nlevel()-1);

		level_threads.emplace_back();
		level_planes.emplace_back();
		level_teams.emplace_back();
		level_map[so.shape(0)] = level_planes.size() - 1;
		this->setup_impl(so, level_planes.back(), level_threads.back(), level_teams.back());
	}
	void run(const stencil_op<seven_pt> & so, grid_func & x,
	         const grid_func & b, cycle::Dir cycle_dir) override { this->run_impl(so, x, b, cycle_dir); }
	void run(const stencil_op<xxvii_pt> & so, grid_func & x,
	         const grid_func & b, cycle::Dir cycle_dir) override { this->run_impl(so, x, b, cycle_dir); }

	template<class sten3, class sten2>
	void setup_impl(const stencil_op<sten3> & so, std::vector<slv2_ptr<sten2>> & planes,
	                std::array<plane_ult<sten2>, 2> & threads,
	                std::array<plane_team, 2> & teams);

	template<class sten>
	void run_impl(const stencil_op<sten> & so, grid_func & x,
	              const grid_func & b, cycle::Dir cycle_dir)
	{
		auto & halo_service = this->services->template get<halo_exchange>();
		if (std::is_same<sten, seven_pt>::value) {
			#ifdef PLANE_AGG
			if (aggregate)
				relax_planes_agg<rdir>(so, x, b, cycle_dir, halo_service, fine_planes, fine_threads);
			else
				relax_planes<rdir>(so, x, b, cycle_dir, halo_service, fine_planes);
			#else
			relax_planes<rdir>(so, x, b, cycle_dir, halo_service, fine_planes);
			#endif
		} else {
			len_t lsize = so.shape(0);
			#ifdef PLANE_AGG
			if (aggregate)
				relax_planes_agg<rdir>(so, x, b, cycle_dir, halo_service, level_planes[level_map[lsize]], level_threads[level_map[lsize]]);
			else
				relax_planes<rdir>(so, x, b, cycle_dir, halo_service, level_planes[level_map[lsize]]);
			#else
			relax_planes<rdir>(so, x, b, cycle_dir, halo_service, level_planes[level_map[lsize]]);
			#endif
		}
	}

protected:
	std::vector<std::vector<slv2_ptr<cdr2::nine_pt>>> level_planes;
	std::vector<slv2_ptr<cdr2::five_pt>> fine_planes;
	std::vector<std::array<plane_ult<cdr2::nine_pt>, 2>> level_threads;
	std::array<plane_ult<cdr2::five_pt>, 2> fine_threads;
	std::map<len_t, std::size_t> level_map;
	std::vector<std::array<plane_team, 2>> level_teams;
	bool aggregate;
	std::shared_ptr<grid_topo> slice_topo(const grid_topo & topo3);
};

}}}

#endif
