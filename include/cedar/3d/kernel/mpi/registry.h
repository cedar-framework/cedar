#ifndef CEDAR_3D_KERNEL_MPI_REGISTRY_H
#define CEDAR_3D_KERNEL_MPI_REGISTRY_H

#include <cedar/types.h>
#include <cedar/mpi_registry.h>
#include <cedar/cycle/types.h>

#include <cedar/3d/mpi/types.h>

#include <cedar/3d/matvec.h>
#include <cedar/3d/mpi/halo.h>
#include <cedar/3d/relax_stencil.h>
#include <cedar/3d/mpi/grid_func.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/3d/mpi/grid_func.h>
#include <cedar/3d/solver.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/residual.h>
#include <cedar/3d/matvec.h>
#include <cedar/3d/inter/mpi/setup_interp.h>
#include <cedar/3d/inter/mpi/interp.h>
#include <cedar/3d/inter/mpi/restrict.h>
#include <cedar/3d/inter/mpi/galerkin_prod.h>
#include <cedar/3d/relax/mpi/setup_relax.h>
#include <cedar/3d/relax/mpi/relax.h>
#include <cedar/3d/cg/setup_cg_redist.h>
#include <cedar/3d/cg/mpi/setup_cg_lu.h>
#include <cedar/3d/cg/mpi/solve_cg.h>
#include <cedar/3d/kernel/setup_nog.h>

namespace cedar { namespace cdr3 { namespace mpi {
			class redist_solver;
}}}


namespace cedar { namespace cdr3 { namespace kernel { namespace mpi {
namespace mpi = cedar::cdr3::mpi;

class registry : public mpi_registry<registry, cdr3::mpi::stypes, cdr3::mpi::redist_solver, solver<xxvii_pt>>
{
public:
	using parent = mpi_registry<registry, cdr3::mpi::stypes, cdr3::mpi::redist_solver, solver<xxvii_pt>>;
registry(std::shared_ptr<kernel_params> params): parent::mpi_registry(params) {}
registry(config::reader & conf) : parent::mpi_registry(conf) {}

	using parent::halo_exchange;

	void setup_nog(grid_topo & topo,
	               len_t min_coarse, int *nog)
	{
		impls::fortran_setup_nog(*params, topo, min_coarse, nog);
	}


	template <class sten>
		void setup_interp(const mpi::stencil_op<sten> & fop,
		                  const mpi::stencil_op<xxvii_pt> & cop,
		                  inter::mpi::prolong_op & P)
	{
		impls::mpi_setup_interp(*params, halof.get(), fop, cop, P);
	}


	template <class sten>
		void galerkin_prod(const inter::mpi::prolong_op & P,
		                   const mpi::stencil_op<sten> & fop,
		                   mpi::stencil_op<xxvii_pt> & cop)
	{
		impls::mpi_galerkin_prod(*params, halof.get(), P, fop, cop);
	}


	template <class sten>
		void setup_relax(const mpi::stencil_op<sten> & so,
		                 relax_stencil & sor)
	{
		impls::mpi_setup_rbgs_point(*params, so, sor);
	}


	template<class sten>
		void relax(const mpi::stencil_op<sten> & so,
		           mpi::grid_func & x,
		           const mpi::grid_func & b,
		           const relax_stencil & sor,
		           cycle::Dir cdir)
	{
		impls::mpi_relax_rbgs_point(*params, halof.get(), so, x, b, sor, cdir);
	}


	void interp_add(const inter::mpi::prolong_op & P,
	                const mpi::grid_func & coarse,
	                const mpi::grid_func & residual,
	                mpi::grid_func & fine)
	{
		impls::mpi_fortran_interp(*params, halof.get(), P, coarse, residual, fine);
	}


	void matvec(const inter::mpi::restrict_op & R,
	            const mpi::grid_func & x,
	            mpi::grid_func & y)
	{
		impls::mpi_fortran_restrict(*params, R, x, y);
	}


	template <class stencil>
	void residual(const mpi::stencil_op<stencil> & so,
	              const mpi::grid_func & x,
	              const mpi::grid_func & b,
	              mpi::grid_func & r)
	{
		impls::mpi_residual_fortran(*params, halof.get(), so, x, b, r);
	}


	template<class sten>
	void matvec(const mpi::stencil_op<sten> & so,
	            const mpi::grid_func & x,
	            mpi::grid_func & y)
	{
		impls::matvec(*params, halof.get(), so, x, y);
	}


	std::unique_ptr<halo_exchanger> halo_create(grid_topo & topo)
	{
		return impls::setup_msg(*params, topo);
	}


	void halo_exchange(mpi::grid_func & f)
	{
		impls::msg_exchange(*params, halof.get(), f);
	}


	template <class sten>
		void halo_stencil_exchange(mpi::stencil_op<sten> & so)
	{
		impls::msg_stencil_exchange(*params, halof.get(), so);
	}


	template<class sten>
		void setup_cg_lu(const mpi::stencil_op<sten> & so,
		                 mpi::grid_func & ABD)
	{
		impls::mpi_setup_cg_lu(*params, halof.get(), so, ABD);
	}


	void solve_cg(mpi::grid_func &x,
	              const mpi::grid_func &b,
	              const mpi::grid_func &ABD,
	              real_t *bbd)
	{
		impls::mpi_solve_cg_lu(*params, halof.get(), x, b, ABD, bbd);
	}


	template <class sten>
	void setup_cg_redist(const mpi::stencil_op<sten> & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<mpi::redist_solver> * bmg,
	                     std::vector<int> & nblocks)
	{
		impls::setup_cg_redist(*params, halof.get(), so, conf, bmg, nblocks);
	}


	void solve_cg_redist(const mpi::redist_solver &bmg,
	                     mpi::grid_func & x,
	                     const mpi::grid_func & b)
	{
		impls::solve_cg_redist(*params, bmg, x, b);
	}
};

}}}}

#endif
