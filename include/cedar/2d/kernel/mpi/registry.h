#ifndef CEDAR_2D_KERNEL_MPI_REGISTRY_H
#define CEDAR_2D_KERNEL_MPI_REGISTRY_H


#include <cedar/mpi_registry>
#include <cedar/types.h>
#include <cedar/cycle/types.h>

#include <cedar/2d/relax_stencil.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/solver.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/inter/mpi/prolong_op.h>
#include <cedar/2d/inter/mpi/restrict_op.h>
#include <cedar/2d/mpi/redist_solver.h>

namespace cedar { namespace cdr2 { namespace kernel { namespace mpi {
namespace mpi = cedar::cdr2::mpi;
class registry : public mpi_registry<mpi::registry, mpi::stencil_op, relax_stencil,
	inter::mpi::prolong_op, inter::mpi::restrict_op, mpi::grid_func, mpi::redist_solver, solver<nine_pt>>
{
public:
	using parent = mpi_registry<mpi::registry, mpi::stencil_op, relax_stencil,
		inter::mpi::prolong_op, inter::mpi::restrict_op, mpi::grid_func, mpi::redist_solver, solver<nine_pt>>;
registry(std::shared_ptr<kernel_params> params): parent::mpi_registry(params) {}
registry(config::reader & conf) : parent::mpi_registry(conf) {}

	template <class sten>
		void setup_interp(const stencil_op<sten> & fop,
		                  const stencil_op<nine_pt> & cop,
		                  inter::prolong_op & P)
	{
		impls::mpi_setup_interp(*params, fop, cop, P);
	}

	template <class sten>
		void galerkin_prod(const inter::prolong_op & P,
		                   const stencil_op<sten> & fop,
		                   stencil_op<nine_pt> & cop)
	{
		impls::mpi_galerkin_prod(*params, P, fop, cop);
	}

	template <class sten>
		void setup_relax(const stencil_op<sten> & so,
		                 relax_stencil & sor)
	{
		impls::mpi_setup_rbgs_point(*params, so, sor);
	}

	template <class sten>
		void setup_relax_x(const stencil_op<sten> & so,
		                   relax_stencil & sor)
	{
		impls::mpi_setup_rbgs_x(*params, so, sor);
	}


	template <class sten>
		void setup_relax_y(const stencil_op<sten> & so,
		                   relax_stencil & sor)
	{
		impls::mpi_setup_rbgs_y(*params, so, sor);
	}


	template<class sten>
		void setup_relax_xy(const stencil_op<sten> & so,
		                   relax_stencil & sor)
	{
		// TODO: add this
		//impls::setup_relax_xy(*params, so, sor);
	}


	template<class sten>
		void relax(const stencil_op<sten> & so,
		           grid_func & x,
		           const grid_func & b,
		           const relax_stencil & sor,
		           cycle::Dir cdir)
	{
		impls::mpi_relax_rbgs_point(*params, so, x, b, sor, cdir);
	}

	template<class sten>
		void setup_cg_lu(const stencil_op<sten> & so,
		                 grid_func & ABD)
	{
		impls::mpi_setup_cg_lu(*params, so, ABD);
	}


	template<class sten>
	void relax_lines_x(const stencil_op<sten> & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func &res,
	                   cycle::Dir cdir)
	{
		impls::mpi_relax_lines_x(*params, so, x, b, sor, res, cdir);
	}


	template<class sten>
		void relax_lines_y(const stencil_op<sten> & so,
		                   grid_func & x,
		                   const grid_func & b,
		                   const relax_stencil & sor,
		                   grid_func &res,
		                   cycle::Dir cdir)
	{
		impls::mpi_relax_lines_y(*params, so, x, b, sor, res, cdir);
	}


	void solve_cg(grid_func &x,
	              const grid_func &b,
	              const grid_func &ABD,
	              real_t *bbd)
	{
		impls::mpi_solve_cg_lu(*params, x, b, ABD, bbd);
	}


	template<class sten>
	void matvec(const stencil_op<sten> & so,
	            const grid_func & x,
	            grid_func & y)
	{
		impls::matvec(so, x, y);
	}


	void matvec(const inter::restrict_op & R,
	            const grid_func & x,
	            grid_func & y)
	{
		impls::mpi_fortran_restrict(*params, R, x, y);
	}


	template <class stencil>
	void residual(const stencil_op<stencil> & so,
	              const grid_func & x,
	              const grid_func & b,
	              grid_func & r)
	{
		impls::mpi_residual_fortran(*params, so, x, b, r);
	}


	void interp_add(const inter::prolong_op & P,
	                const grid_func & coarse,
	                const grid_func & residual,
	                grid_func & fine)
	{
		impls::mpi_fortran_interp(*params, P, coarse, residual, fine);
	}


	void halo_setup(grid_topo & topo,
	                void **halo_ctx)
	{
		impls::setup_msg(*params, topo, halo_ctx);
	}


	void halo_exchange(grid_func & f)
	{
		impls::msg_exchange(*params, f);
	}


	template <class sten>
		void halo_stencil_exchange(stencil_op<sten> & so)
	{
		impls::msg_stencil_exchange(*params, so);
	}


	template <class sten>
	void setup_cg_boxmg(const stencil_op<sten> & so,
	                    std::shared_ptr<config::reader> conf,
	                    std::shared_ptr<serial_solver> *bmg)
	{
		impls::setup_cg_boxmg(*params, so, conf, bmg);
	}

	void solve_cg_boxmg(const serial_solver &bmg,
	                    grid_func &x,
	                    const grid_func &b)
	{
		impls::solve_cg_boxmg(*params, bmg, x, b);
	}


	template <class sten>
	void setup_cg_redist(const stencil_op<sten> & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<redist_solver> * bmg,
	                     std::vector<int> & nblocks)
	{
		impls::setup_cg_redist(*params, so, conf, bmg, nblocks);
	}


	void solve_cg_redist(const redist_solver &bmg,
	                     grid_func & x,
	                     const grid_func & b)
	{
		impls::solve_cg_redist(*params, bmg, x, b);
	}
};

}}}}


#endif
