#ifndef BOXMG_MPI_REGISTRY_H
#define BOXMG_MPI_REGISTRY_H

#include <boxmg/kernel_registry.h>
#include <boxmg/mpi/grid_topo.h>
#include <boxmg/config/reader.h>


namespace boxmg {

template <class stencil_op, class relax_stencil, class prolong_op, class grid_func, class redist_solver, class serial_solver>

class mpi_registry : public kernel_registry<stencil_op, relax_stencil, prolong_op, grid_func>
{
public:
	using kernel_registry<stencil_op, relax_stencil, prolong_op, grid_func>::active;
	void setup_nog(grid_topo &topo,
	               len_t min_coarse, int *nog)
	{
		active.run(kernel_name::setup_nog, topo,
		           static_cast<len_t>(min_coarse),
		           static_cast<int*>(nog));
	}


	void halo_setup(grid_topo &topo,
	                void **halo_ctx)
	{
		active.run(kernel_name::halo_setup,
		           static_cast<grid_topo&>(topo),
		           static_cast<void**>(halo_ctx));
	}


	void halo_exchange(grid_func &f)
	{
		active.run(kernel_name::halo_exchange,
		           static_cast<grid_func&>(f));
	}


	void halo_exchange(const grid_func &f, void *halo_ctx)
	{
		auto & fd = const_cast<grid_func&>(f);
		fd.halo_ctx = halo_ctx;
		active.run(kernel_name::halo_exchange,
		           static_cast<grid_func&>(fd));

	}


	void halo_stencil_exchange(stencil_op & so)
	{
		active.run(kernel_name::halo_stencil_exchange,
		           so);
	}


	void setup_cg_boxmg(const stencil_op & so,
	                    std::shared_ptr<config::reader> conf,
	                    std::shared_ptr<serial_solver> *bmg)
	{
		active.run(kernel_name::setup_cg_boxmg, so,
		           static_cast<std::shared_ptr<config::reader>>(conf),
		           static_cast<std::shared_ptr<serial_solver>*>(bmg));

	}


	void solve_cg_boxmg(const serial_solver &bmg,
	                    grid_func &x,
	                    const grid_func &b)
	{
		active.run(kernel_name::solve_cg_boxmg, bmg, x, b);
	}


	void setup_cg_redist(const stencil_op & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<redist_solver> * bmg,
	                     std::vector<int> & nblocks)
	{
		active.run(kernel_name::setup_cg_redist, so,
		           static_cast<std::shared_ptr<config::reader>>(conf),
		           static_cast<std::shared_ptr<redist_solver>*>(bmg), nblocks);
	}


	void solve_cg_redist(const redist_solver &bmg,
	                     grid_func & x,
	                     const grid_func & b)
	{
		active.run(kernel_name::solve_cg_redist,
		           bmg, x, b);
	}


	void matvec(const stencil_op & so,
	            const grid_func &x,
	            grid_func &b)
	{
		active.run(kernel_name::matvec, so, x, b);
	}
};

}
#endif
