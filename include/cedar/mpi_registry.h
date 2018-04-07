#ifndef CEDAR_MPI_REGISTRY_H
#define CEDAR_MPI_REGISTRY_H

#include <cedar/kernel_registry.h>
#include <cedar/mpi/grid_topo.h>

namespace cedar {


	/**
	   Base class used to define kernels used in a distributed
	   multilevel solve.

	   By subclassing this, you must implement each method to provide
	   a registry of implementations to the common distributed
	   multilevel operations.  This class given to multilevel solver
	   objects so these operations can be performed generically.

	   @tparam child The class that is inheriting this base.
	   @tparam solver_types A structured listing of types used in a multilevel solve.
	   @tparam redist_solver Class for a redistributed coarse-grid solver.
	   @tparam serial_solver Class for a serial coarse-grid solve.
	*/
	template <class child,
		class solver_types,
		class redist_solver,
		class serial_solver,
		class halo_exchanger>
class mpi_registry : public kernel_registry<child,
                                    		solver_types>
{
public:
	using parent = kernel_registry<child, solver_types>;
	// extract basic types from solver_types
	template<class sten>
	using stencil_op = typename solver_types::template stencil_op<sten>;
	using grid_func = typename solver_types::grid_func;
	using prolong_op = typename solver_types::prolong_op;
	using restrict_op = typename solver_types::restrict_op;
	using relax_stencil = typename solver_types::relax_stencil;

mpi_registry(std::shared_ptr<kernel_params> params): parent::kernel_registry(params) {
	}
mpi_registry(config::reader & conf) : parent::kernel_registry(conf) {
	}

	void setup_nog(grid_topo &topo,
	               len_t min_coarse, int *nog)
	{
		log::debug << "Running kernel <setup_nog>" << std::endl;
		static_cast<child*>(this)->setup_nog(topo, min_coarse, nog);
	}


	void halo_setup(grid_topo &topo)
	{
		log::debug << "Running kernel <halo_setup>" << std::endl;
		//halof = static_cast<child*>(this)->halo_create(topo);
		halof = std::make_shared<halo_exchanger>(*(this->params), topo);
	}


	void halo_exchange(grid_func &f)
	{
		log::debug << "Running kernel <halo_exchange>" << std::endl;
		//static_cast<child*>(this)->halo_exchange(f);
		halof->exchange(f);
	}


	template <class sten>
	void halo_stencil_exchange(stencil_op<sten> & so)
	{
		log::debug << "Running kernel <halo_stencil_exchange>" << std::endl;
		//static_cast<child*>(this)->halo_stencil_exchange(so);
		halof->exchange(so);
	}


	template <class sten>
	void setup_cg_boxmg(const stencil_op<sten> & so,
	                    std::shared_ptr<config::reader> conf,
	                    std::shared_ptr<serial_solver> *bmg)
	{
		log::debug << "Running kernel <setup_cg_boxmg>" << std::endl;
		static_cast<child*>(this)->setup_cg_boxmg(so, conf, bmg);
	}


	void solve_cg_boxmg(const serial_solver &bmg,
	                    grid_func &x,
	                    const grid_func &b)
	{
		log::debug << "Running kernel <solve_cg_boxmg>" << std::endl;
		static_cast<child*>(this)->solve_cg_boxmg(bmg, x, b);
	}


	template <class sten>
	void setup_cg_redist(const stencil_op<sten> & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<redist_solver> * bmg,
	                     std::vector<int> & nblocks)
	{
		log::debug << "Running kernel <setup_cg_redist>" << std::endl;
		static_cast<child*>(this)->setup_cg_redist(so, conf, bmg, nblocks);
	}


	void solve_cg_redist(const redist_solver &bmg,
	                     grid_func & x,
	                     const grid_func & b)
	{
		log::debug << "Running kernel <solve_cg_redist>" << std::endl;
		static_cast<child*>(this)->solve_cg_redist(bmg, x, b);
	}

	std::shared_ptr<halo_exchanger> get_halo_exchanger() { return halof; }

protected:
	std::shared_ptr<halo_exchanger> halof;
};

}
#endif
