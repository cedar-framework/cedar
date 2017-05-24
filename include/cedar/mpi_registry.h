#ifndef CEDAR_MPI_REGISTRY_H
#define CEDAR_MPI_REGISTRY_H

#include <cedar/kernel_registry.h>
#include <cedar/mpi/grid_topo.h>

namespace cedar {

	template <class child,
		template<class> class stencil_op,
		class relax_stencil,
		class prolong_op,
		class restrict_op,
		class grid_func,
		class redist_solver,
		class serial_solver>
class mpi_registry : public kernel_registry<child,
		                                    stencil_op,
		                                    relax_stencil,
		                                    prolong_op,
		                                    restrict_op,
		                                    grid_func>
{
public:
	using parent = kernel_registry<child, stencil_op, relax_stencil, prolong_op, restrict_op, grid_func>;
mpi_registry(std::shared_ptr<kernel_params> params): parent::kernel_registry)(params) {}
mpi_registry(config::reader & conf) : parent::kernel_registry(conf) {}

	void setup_nog(grid_topo &topo,
	               len_t min_coarse, int *nog)
	{
		static_cast<child*>(this)->setup_nog(topo, min_coarse, nog);
	}


	void halo_setup(grid_topo &topo,
	                void **halo_ctx)
	{
		static_cast<child*>(this)->halo_setup(topo, halo_ctx);
	}


	void halo_exchange(grid_func &f)
	{
		static_cast<child*>(this)->halo_exchange(f);
	}


	void halo_exchange(grid_func &f, void *halo_ctx)
	{
		f.halo_ctx = halo_ctx;
		static_cast<child*>(this)->halo_exchange(f);
	}


	template <class sten>
	void halo_stencil_exchange(stencil_op<sten> & so)
	{
		static_cast<child*>(this)->halo_stencil_exchange(so);
	}


	template <class sten>
	void setup_cg_boxmg(const stencil_op<sten> & so,
	                    std::shared_ptr<config::reader> conf,
	                    std::shared_ptr<serial_solver> *bmg)
	{
		static_cast<child*>(this)->setup_cg_boxmg(so, conf, bmg);
	}


	void solve_cg_boxmg(const serial_solver &bmg,
	                    grid_func &x,
	                    const grid_func &b)
	{
		static_cast<child*>(this)->solve_cg_boxmg(bmg, x, b);
	}


	template <class sten>
	void setup_cg_redist(const stencil_op<sten> & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<redist_solver> * bmg,
	                     std::vector<int> & nblocks)
	{
		static_cast<child*>(this)->setup_cg_redist(so, conf, bmg, nblocks);
	}


	void solve_cg_redist(const redist_solver &bmg,
	                     grid_func & x,
	                     const grid_func & b)
	{
		static_cast<child*>(this)->solve_cg_redist(bmg, x, b);
	}
};

}
#endif
