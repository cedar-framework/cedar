#ifndef CEDAR_2D_KERNEL_MPI_REGISTRY_H
#define CEDAR_2D_KERNEL_MPI_REGISTRY_H


#include <cedar/mpi_registry.h>
#include <cedar/types.h>
#include <cedar/cycle/types.h>
#include <cedar/2d/mpi/types.h>

#include <cedar/2d/matvec.h>
#include <cedar/2d/relax/mpi/relax.h>
#include <cedar/2d/relax/mpi/ml_relax.h>
#include <cedar/2d/relax_stencil.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/solver.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/msg_exchanger.h>
#include <cedar/2d/inter/mpi/prolong_op.h>
#include <cedar/2d/inter/mpi/restrict_op.h>
#include <cedar/2d/inter/mpi/galerkin_prod.h>
#include <cedar/2d/cg/setup_cg_boxmg.h>
#include <cedar/2d/cg/mpi/setup_cg_lu.h>
#include <cedar/2d/kernel/setup_nog.h>

namespace cedar { namespace cdr2 { namespace kernel { namespace mpi {
	namespace mpi = cedar::cdr2::mpi;
	template<class halo = mpi::msg_exchanger>
	class registry : public mpi_registry<registry<halo>, cdr2::mpi::stypes, solver<nine_pt>,
		halo>
{
public:
	using parent = mpi_registry<registry<halo>, cdr2::mpi::stypes, solver<nine_pt>, halo>;
registry(std::shared_ptr<kernel_params> params): parent::mpi_registry(params),
     xrelax(mpi::relax_dir::x), yrelax(mpi::relax_dir::y) {}
registry(config::reader & conf) : parent::mpi_registry(conf),
     xrelax(mpi::relax_dir::x), yrelax(mpi::relax_dir::y) {}

	using parent::halo_exchange;
	using parent::halo_setup;
	using parent::halo_stencil_exchange;

	void setup_nog(grid_topo & topo,
	               len_t min_coarse, int *nog)
	{
		impls::fortran_setup_nog(*this->params, topo, min_coarse, nog);
	}

	template <class sten>
		void setup_interp(const mpi::stencil_op<sten> & fop,
		                  const mpi::stencil_op<nine_pt> & cop,
		                  inter::mpi::prolong_op & P)
	{
		impls::mpi_setup_interp(*this->params, this->halof.get(), fop, cop, P);
	}

	template <class sten>
		void galerkin_prod(const inter::mpi::prolong_op & P,
		                   const mpi::stencil_op<sten> & fop,
		                   mpi::stencil_op<nine_pt> & cop)
	{
		impls::mpi_galerkin_prod(*this->params, this->halof.get(), P, fop, cop);
	}

	template <class sten>
		void setup_relax(const mpi::stencil_op<sten> & so,
		                 relax_stencil & sor)
	{
		impls::mpi_setup_rbgs_point(*this->params, so, sor);
	}

	template <class sten>
		void setup_relax_x(const mpi::stencil_op<sten> & so,
		                   relax_stencil & sor)
	{
		if (this->params->ml_relax.enabled) {
			xrelax.setup(*this->params, so, sor);
		} else {
			impls::mpi_setup_rbgs_x(*this->params, so, sor);
		}
	}


	template <class sten>
		void setup_relax_y(const mpi::stencil_op<sten> & so,
		                   relax_stencil & sor)
	{
		if (this->params->ml_relax.enabled) {
			yrelax.setup(*this->params, so, sor);
		} else {
			impls::mpi_setup_rbgs_y(*this->params, so, sor);
		}
	}


	template<class sten>
		void setup_relax_xy(const mpi::stencil_op<sten> & so,
		                   relax_stencil & sor)
	{
		// TODO: add this
		//impls::setup_relax_xy(*this->params, so, sor);
	}


	template<class sten>
		void relax(const mpi::stencil_op<sten> & so,
		           mpi::grid_func & x,
		           const mpi::grid_func & b,
		           const relax_stencil & sor,
		           cycle::Dir cdir)
	{
		impls::mpi_relax_rbgs_point(*this->params, this->halof.get(), so, x, b, sor, cdir);
	}

	template<class sten>
		void setup_cg_lu(const mpi::stencil_op<sten> & so,
		                 mpi::grid_func & ABD)
	{
		log::error << "This kernel (solve_cg) is unavailable for mpi solver!" << std::endl;
		/* impls::mpi_setup_cg_lu(*this->params, this->halof.get(), so, ABD); */
	}


	template<class sten>
	void relax_lines_x(const mpi::stencil_op<sten> & so,
	                   mpi::grid_func & x,
	                   const mpi::grid_func & b,
	                   const relax_stencil & sor,
	                   mpi::grid_func &res,
	                   cycle::Dir cdir)
	{
		if (this->params->ml_relax.enabled) {
			xrelax.solve(*this->params, this->halof.get(), so, x, b, sor, res, cdir);
		} else {
			impls::mpi_relax_lines_x(*this->params, this->halof.get(), so, x, b, sor, res, cdir);
		}
	}


	template<class sten>
		void relax_lines_y(const mpi::stencil_op<sten> & so,
		                   mpi::grid_func & x,
		                   const mpi::grid_func & b,
		                   const relax_stencil & sor,
		                   mpi::grid_func &res,
		                   cycle::Dir cdir)
	{
		if (this->params->ml_relax.enabled) {
			yrelax.solve(*this->params, this->halof.get(), so, x, b, sor, res, cdir);
		} else {
			impls::mpi_relax_lines_y(*this->params, this->halof.get(), so, x, b, sor, res, cdir);
		}
	}


	void solve_cg(mpi::grid_func &x,
	              const mpi::grid_func &b,
	              const mpi::grid_func &ABD,
	              real_t *bbd)
	{
		log::error << "This kernel (solve_cg) is deprecated for mpi solver!" << std::endl;
		/* impls::mpi_solve_cg_lu(*this->params, this->halof.get(), x, b, ABD, bbd); */
	}


	template<class sten>
	void matvec(const mpi::stencil_op<sten> & so,
	            const mpi::grid_func & x,
	            mpi::grid_func & y)
	{
		impls::matvec(*this->params, this->halof.get(), so, x, y);
	}


	void matvec(const inter::mpi::restrict_op & R,
	            const mpi::grid_func & x,
	            mpi::grid_func & y)
	{
		impls::mpi_fortran_restrict(*this->params, R, x, y);
	}


	template <class stencil>
	void residual(const mpi::stencil_op<stencil> & so,
	              const mpi::grid_func & x,
	              const mpi::grid_func & b,
	              mpi::grid_func & r)
	{
		impls::mpi_residual_fortran(*this->params, this->halof.get(), so, x, b, r);
	}


	void interp_add(const inter::mpi::prolong_op & P,
	                const mpi::grid_func & coarse,
	                const mpi::grid_func & residual,
	                mpi::grid_func & fine)
	{
		impls::mpi_fortran_interp(*this->params, this->halof.get(), P, coarse, residual, fine);
	}


	template <class sten>
	void setup_cg_boxmg(const mpi::stencil_op<sten> & so,
	                    std::shared_ptr<config::reader> conf,
	                    std::shared_ptr<solver<nine_pt>> *bmg)
	{
		log::error << "This kernel (setup_cg_boxmg) is deprecated for mpi solver!" << std::endl;
		/* impls::setup_cg_boxmg(*this->params, this->halof.get(), so, conf, bmg); */
	}

	void solve_cg_boxmg(const solver<nine_pt> &bmg,
	                    mpi::grid_func &x,
	                    const mpi::grid_func &b)
	{
		log::error << "This kernel (solve_cg_boxmg) is deprecated for mpi solver!" << std::endl;
		/* impls::solve_cg_boxmg(*this->params, this->halof.get(), bmg, x, b); */
	}

protected:
	mpi::ml_line_relaxer xrelax, yrelax;

};

}}}}


#endif
