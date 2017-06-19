#ifndef CEDAR_MPI_REGISTRY_H
#define CEDAR_MPI_REGISTRY_H

#include <cedar/kernel_registry.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/halo_exchanger.h>

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
		class serial_solver, short ND>
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
		halof.exchange = [this](int k, int nog, real_t *so_data, std::array<len_t, ND> len, void *halo_ctx) {
			this->halo_exchange(k, nog, so_data, len, halo_ctx);
		};
	}
mpi_registry(config::reader & conf) : parent::kernel_registry(conf) {
		halof.exchange = [this](int k, int nog, real_t *so_data, std::array<len_t, ND> len, void *halo_ctx) {
			this->halo_exchange(k, nog, so_data, len, halo_ctx);
		};
	}

	void setup_nog(grid_topo &topo,
	               len_t min_coarse, int *nog)
	{
		log::debug << "Running kernel <setup_nog>" << std::endl;
		static_cast<child*>(this)->setup_nog(topo, min_coarse, nog);
	}


	void halo_setup(grid_topo &topo,
	                void **halo_ctx)
	{
		log::debug << "Running kernel <halo_setup>" << std::endl;
		static_cast<child*>(this)->halo_setup(topo, halo_ctx);
	}


	void halo_exchange(grid_func &f)
	{
		log::debug << "Running kernel <halo_exchange>" << std::endl;
		static_cast<child*>(this)->halo_exchange(f);
	}


	void halo_exchange(grid_func &f, void *halo_ctx)
	{
		log::debug << "Running kernel <halo_exchange>" << std::endl;
		f.halo_ctx = halo_ctx;
		static_cast<child*>(this)->halo_exchange(f);
	}


	void halo_exchange(int k, int nog, real_t* so_data, std::array<len_t, ND> len, void *halo_ctx)
	{
		log::debug << "Running kernel <halo_exchange>" << std::endl;
		static_cast<child*>(this)->halo_exchange(k, nog, so_data, len, halo_ctx);
	}


	template <class sten>
	void halo_stencil_exchange(stencil_op<sten> & so)
	{
		log::debug << "Running kernel <halo_stencil_exchange>" << std::endl;
		static_cast<child*>(this)->halo_stencil_exchange(so);
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

protected:
	halo_exchanger<ND> halof;
};

}
#endif
