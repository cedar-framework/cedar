#ifndef CEDAR_3D_KERNEL_REGISTRY_H
#define CEDAR_3D_KERNEL_REGISTRY_H


#include <cedar/kernel_registry.h>
#include <cedar/types.h>
#include <cedar/3d/types.h>
#include <cedar/cycle/types.h>
#include <cedar/3d/inter/setup_interp.h>
#include <cedar/3d/inter/galerkin_prod.h>
#include <cedar/3d/relax/setup_relax.h>
#include <cedar/3d/cg/setup_cg_lu.h>
#include <cedar/3d/relax/relax.h>
#include <cedar/3d/cg/solve_cg.h>
#include <cedar/3d/inter/restrict.h>
#include <cedar/3d/inter/interp.h>
#include <cedar/3d/residual.h>

namespace cedar { namespace cdr3 { namespace kernel {

class registry : public kernel_registry<registry, stypes>
{
public:
	using parent = kernel_registry<registry, stypes>;
registry(std::shared_ptr<kernel_params> params) : parent::kernel_registry(params) {}
registry(config::reader & conf) : parent::kernel_registry(conf) {}

	template<class sten>
		void setup_interp(const stencil_op<sten> & fop,
		                  const stencil_op<xxvii_pt> & cop,
		                  inter::prolong_op & P)
	{
		impls::setup_interp(*params, fop, cop, P);
	}

	template<class sten>
		void galerkin_prod(const inter::prolong_op & P,
		                   const stencil_op<sten> & fop,
		                   stencil_op<xxvii_pt> & cop)
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
		impls::residual(*params, so, x, b, r);
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
