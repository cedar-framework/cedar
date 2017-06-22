#ifndef CEDAR_2D_KERNEL_MPI_REGISTRY_H
#define CEDAR_2D_KERNEL_MPI_REGISTRY_H


#include <cedar/mpi_registry.h>
#include <cedar/types.h>
#include <cedar/cycle/types.h>
#include <cedar/2d/mpi/types.h>

#include <cedar/2d/matvec.h>
#include <cedar/2d/mpi/halo.h>
#include <cedar/2d/relax/mpi/relax.h>
#include <cedar/2d/relax_stencil.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/solver.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/inter/mpi/prolong_op.h>
#include <cedar/2d/inter/mpi/restrict_op.h>
#include <cedar/2d/inter/mpi/galerkin_prod.h>
#include <cedar/2d/cg/setup_cg_boxmg.h>
#include <cedar/2d/cg/setup_cg_redist.h>
#include <cedar/2d/kernel/setup_nog.h>

namespace cedar { namespace cdr2 { namespace mpi {
			class redist_solver;
		}
	}
}

namespace cedar { namespace cdr2 { namespace kernel { namespace mpi {
	namespace mpi = cedar::cdr2::mpi;
	class registry : public mpi_registry<registry, cdr2::mpi::stypes, cdr2::mpi::redist_solver, solver<nine_pt>, 2>
{
public:
	using parent = mpi_registry<registry, cdr2::mpi::stypes, cdr2::mpi::redist_solver, solver<nine_pt>, 2>;
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
		                  const mpi::stencil_op<nine_pt> & cop,
		                  inter::mpi::prolong_op & P)
	{
		impls::mpi_setup_interp(*params, halof, fop, cop, P);
	}

	template <class sten>
		void galerkin_prod(const inter::mpi::prolong_op & P,
		                   const mpi::stencil_op<sten> & fop,
		                   mpi::stencil_op<nine_pt> & cop)
	{
		impls::mpi_galerkin_prod(*params, halof, P, fop, cop);
	}

	template <class sten>
		void setup_relax(const mpi::stencil_op<sten> & so,
		                 relax_stencil & sor)
	{
		impls::mpi_setup_rbgs_point(*params, so, sor);
	}

	template <class sten>
		void setup_relax_x(const mpi::stencil_op<sten> & so,
		                   relax_stencil & sor)
	{
		impls::mpi_setup_rbgs_x(*params, so, sor);
	}


	template <class sten>
		void setup_relax_y(const mpi::stencil_op<sten> & so,
		                   relax_stencil & sor)
	{
		impls::mpi_setup_rbgs_y(*params, so, sor);
	}


	template<class sten>
		void setup_relax_xy(const mpi::stencil_op<sten> & so,
		                   relax_stencil & sor)
	{
		// TODO: add this
		//impls::setup_relax_xy(*params, so, sor);
	}


	template<class sten>
		void relax(const mpi::stencil_op<sten> & so,
		           mpi::grid_func & x,
		           const mpi::grid_func & b,
		           const relax_stencil & sor,
		           cycle::Dir cdir)
	{
		impls::mpi_relax_rbgs_point(*params, halof, so, x, b, sor, cdir);
	}

	template<class sten>
		void setup_cg_lu(const mpi::stencil_op<sten> & so,
		                 mpi::grid_func & ABD)
	{
		impls::mpi_setup_cg_lu(*params, so, ABD);
	}


	template<class sten>
	void relax_lines_x(const mpi::stencil_op<sten> & so,
	                   mpi::grid_func & x,
	                   const mpi::grid_func & b,
	                   const relax_stencil & sor,
	                   mpi::grid_func &res,
	                   cycle::Dir cdir)
	{
		impls::mpi_relax_lines_x(*params, so, x, b, sor, res, cdir);
	}


	template<class sten>
		void relax_lines_y(const mpi::stencil_op<sten> & so,
		                   mpi::grid_func & x,
		                   const mpi::grid_func & b,
		                   const relax_stencil & sor,
		                   mpi::grid_func &res,
		                   cycle::Dir cdir)
	{
		impls::mpi_relax_lines_y(*params, so, x, b, sor, res, cdir);
	}


	void solve_cg(mpi::grid_func &x,
	              const mpi::grid_func &b,
	              const mpi::grid_func &ABD,
	              real_t *bbd)
	{
		impls::mpi_solve_cg_lu(*params, x, b, ABD, bbd);
	}


	template<class sten>
	void matvec(const mpi::stencil_op<sten> & so,
	            const mpi::grid_func & x,
	            mpi::grid_func & y)
	{
		impls::matvec(*params, so, x, y);
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
		impls::mpi_residual_fortran(*params, halof, so, x, b, r);
	}


	void interp_add(const inter::mpi::prolong_op & P,
	                const mpi::grid_func & coarse,
	                const mpi::grid_func & residual,
	                mpi::grid_func & fine)
	{
		impls::mpi_fortran_interp(*params, P, coarse, residual, fine);
	}


	std::unique_ptr<halo_exchanger> halo_create(grid_topo & topo)
	{
		return impls::setup_msg(*params, topo);
	}


	void halo_exchange(mpi::grid_func & f)
	{
		impls::msg_exchange(*params, f);
	}


	template <class sten>
		void halo_stencil_exchange(mpi::stencil_op<sten> & so)
	{
		impls::msg_stencil_exchange(*params, so);
	}


	template <class sten>
	void setup_cg_boxmg(const mpi::stencil_op<sten> & so,
	                    std::shared_ptr<config::reader> conf,
	                    std::shared_ptr<solver<nine_pt>> *bmg)
	{
		impls::setup_cg_boxmg(*params, so, conf, bmg);
	}

	void solve_cg_boxmg(const solver<nine_pt> &bmg,
	                    mpi::grid_func &x,
	                    const mpi::grid_func &b)
	{
		impls::solve_cg_boxmg(*params, bmg, x, b);
	}


	template <class sten>
	void setup_cg_redist(const mpi::stencil_op<sten> & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<mpi::redist_solver> * bmg,
	                     std::vector<int> & nblocks)
	{
		impls::setup_cg_redist(*params, so, conf, bmg, nblocks);
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
