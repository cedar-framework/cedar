#ifndef CEDAR_2D_KERNEL_REGISTRY_H
#define CEDAR_2D_KERNEL_REGISTRY_H


#include <cedar/kernel_registry.h>
#include <cedar/types.h>
#include <cedar/cycle/types.h>

#include <cedar/2d/relax_stencil.h>
#include <cedar/2d/grid_func.h>
#include <cedar/2d/stencil_op.h>

#include <cedar/2d/inter/restrict_op.h>
#include <cedar/2d/inter/prolong_op.h>
#include <cedar/2d/inter/restrict.h>
#include <cedar/2d/inter/interp.h>
#include <cedar/2d/inter/setup_interp.h>
#include <cedar/2d/inter/galerkin_prod.h>
#include <cedar/2d/relax/setup_relax.h>
#include <cedar/2d/relax/relax.h>
#include <cedar/2d/cg/setup_cg_lu.h>
#include <cedar/2d/cg/solve_cg.h>
#include <cedar/2d/residual.h>


namespace cedar { namespace cdr2 { namespace kernel {

			class registry : public kernel_registry<registry, stencil_op, relax_stencil,
				inter::prolong_op, inter::restrict_op, grid_func>
{
public:
	using parent = kernel_registry<registry, stencil_op, relax_stencil,
		inter::prolong_op, inter::restrict_op, grid_func>;
registry(std::shared_ptr<kernel_params> params) : parent::kernel_registry(params) {}
registry(config::reader & conf) : parent::kernel_registry(conf) {}
	template<class sten>
		void setup_interp(const stencil_op<sten> & fop,
		                  const stencil_op<nine_pt> & cop,
		                  inter::prolong_op & P)
	{
		impls::setup_interp(*params, fop, cop, P);
	}

	template<class sten>
		void galerkin_prod(const inter::prolong_op & P,
		                   const stencil_op<sten> & fop,
		                   stencil_op<nine_pt> & cop)
	{
		impls::galerkin_prod(*params, P, fop, cop);
	}

	template<class sten>
		void setup_relax(const stencil_op<sten> & so,
		                 relax_stencil & sor)
	{
		impls::setup_rbgs_point(*params, so, sor);
	}

	template<class sten>
		void setup_relax_x(const stencil_op<sten> & so,
		                   relax_stencil & sor)
	{
		impls::setup_rbgs_x(*params, so, sor);
	}

	template<class sten>
		void setup_relax_y(const stencil_op<sten> & so,
		                   relax_stencil & sor)
	{
		impls::setup_rbgs_y(*params, so, sor);
	}

	template<class sten>
		void setup_relax_xy(const stencil_op<sten> & so,
		                   relax_stencil & sor)
	{
		// TODO: add this
		//impls::setup_relax_xy(*params, so, sor);
	}

	template<class sten>
		void setup_cg_lu(const stencil_op<sten> & so,
		                 grid_func & ABD)
	{
		impls::setup_cg_lu(*params, so, ABD);
	}

	template<class sten>
		void relax(const stencil_op<sten> & so,
		           grid_func & x,
		           const grid_func & b,
		           const relax_stencil & sor,
		           cycle::Dir cdir)
	{
		impls::relax_rbgs_point(*params, so, x, b, sor, cdir);
	}

	template<class sten>
	void relax_lines_x(const stencil_op<sten> & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func &res,
	                   cycle::Dir cdir)
	{
		impls::relax_lines_x(*params, so, x, b, sor, res, cdir);
	}

	template<class sten>
		void relax_lines_y(const stencil_op<sten> & so,
		                   grid_func & x,
		                   const grid_func & b,
		                   const relax_stencil & sor,
		                   grid_func &res,
		                   cycle::Dir cdir)
	{
		impls::relax_lines_y(*params, so, x, b, sor, res, cdir);
	}

	void solve_cg(grid_func &x,
	              const grid_func &b,
	              const grid_func &ABD,
	              real_t *bbd)
	{
		impls::fortran_solve_cg(*params, x, b, ABD, bbd);
	}


	template<class sten>
	void matvec(const stencil_op<sten> & so,
	            const grid_func & x,
	            grid_func & y)
	{
		log::error << "matvec not defined!" << std::endl;
	}


	void matvec(const inter::restrict_op & R,
	            const grid_func & x,
	            grid_func & y)
	{
		impls::fortran_restrict(*params, R, x, y);
	}


	template <class stencil>
	void residual(const stencil_op<stencil> & so,
	              const grid_func & x,
	              const grid_func & b,
	              grid_func & r)
	{
		impls::residual_fortran(*params, so, x, b, r);
	}


	void interp_add(const inter::prolong_op & P,
	                const grid_func & coarse,
	                const grid_func & residual,
	                grid_func & fine)
	{
		impls::fortran_interp(*params, P, coarse, residual, fine);
	}

};

}}}


#endif
