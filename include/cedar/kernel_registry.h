#ifndef CEDAR_KERNEL_REGISTRY_H
#define CEDAR_KERNEL_REGISTRY_H

#include <memory>

#include <cedar/cycle/types.h>
#include <cedar/types.h>
#include <cedar/kernel_params.h>
#include <cedar/config/reader.h>

namespace cedar {


	template <class child,
		class solver_types>
struct kernel_registry
{
	// extract basic types from solver_types
	template<class sten>
	using stencil_op = typename solver_types::template stencil_op<sten>;
	using grid_func = typename solver_types::grid_func;
	using prolong_op = typename solver_types::prolong_op;
	using restrict_op = typename solver_types::restrict_op;
	using relax_stencil = typename solver_types::relax_stencil;

    kernel_registry(std::shared_ptr<kernel_params> params) : params(params) {}
	kernel_registry(config::reader & conf) { params = build_kernel_params(conf); }

	template<class stencil0, class stencil1>
	void setup_interp(const stencil_op<stencil0> & fop,
	                  const stencil_op<stencil1> &cop, prolong_op & P) {
		log::debug << "Running kernel <setup_interp>" << std::endl;
		static_cast<child*>(this)->setup_interp(fop, cop, P);
	}


	template<class stencil0, class stencil1>
	void galerkin_prod(const prolong_op & P,
	                   const stencil_op<stencil0> & fop,
	                   stencil_op<stencil1> & cop) {
		log::debug << "Running kernel <galerkin_prod>" << std::endl;
		static_cast<child*>(this)->galerkin_prod(P, fop, cop);
	}


	template<class stencil>
	void setup_relax(const stencil_op<stencil> & so,
	                 relax_stencil & sor) {
		log::debug << "Running kernel <setup_relax>" << std::endl;
		static_cast<child*>(this)->setup_relax(so, sor);
	}


	template<class stencil>
	void setup_relax_x(const stencil_op<stencil> & so,
	                   relax_stencil & sor) {
		log::debug << "Running kernel <setup_relax_x>" << std::endl;
		static_cast<child*>(this)->setup_relax_x(so, sor);
	}


	template<class stencil>
	void setup_relax_y(const stencil_op<stencil> & so,
	                   relax_stencil & sor) {
		log::debug << "Running kernel <setup_relax_y>" << std::endl;
		static_cast<child*>(this)->setup_relax_y(so, sor);
	}


	template<class stencil>
	void setup_relax_xy(const stencil_op<stencil> & so,
	                            relax_stencil & sor) {
		log::debug << "Running kernel <setup_relax_xy" << std::endl;
		static_cast<child*>(this)->setup_relax_xy(so, sor);
	}


	template<class stencil>
	void setup_cg_lu(const stencil_op<stencil> & so,
	                 grid_func & ABD) {
		log::debug << "Running kernel <setup_cg_lu>" << std::endl;
		static_cast<child*>(this)->setup_cg_lu(so, ABD);
	}


	template<class stencil>
	void relax(const stencil_op<stencil> & so,
	           grid_func & x,
	           const grid_func & b,
	           const relax_stencil & sor,
	           cycle::Dir cdir) {
		log::debug << "Running kernel <relax>" << std::endl;
		static_cast<child*>(this)->relax(so, x, b, sor, cdir);
	}


	template<class stencil>
	void relax_lines_x(const stencil_op<stencil> & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func &res,
	                   cycle::Dir cdir) {
		log::debug << "Running kernel <relax_lines_x>" << std::endl;
		static_cast<child*>(this)->relax_lines_x(so, x, b, sor, res, cdir);
	}


	template<class stencil>
	void relax_lines_y(const stencil_op<stencil> & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func &res,
	                   cycle::Dir cdir) {
		log::debug << "Running kernel <relax_lines_y>" << std::endl;
		static_cast<child*>(this)->relax_lines_y(so, x, b, sor, res, cdir);
	}


	void solve_cg(grid_func &x,
	              const grid_func &b,
	              const grid_func &ABD,
	              real_t *bbd) {
		log::debug << "Running kernel <solve_cg>" << std::endl;
		static_cast<child*>(this)->solve_cg(x, b, ABD, bbd);
	}


	template<class stencil>
	void matvec(const stencil_op<stencil> & so,
	            const grid_func & x,
	            grid_func & y) {
		log::debug << "Running kernel <matvec>" << std::endl;
		static_cast<child*>(this)->matvec(so, x, y);
	}


	void matvec(const restrict_op & R,
	            const grid_func & x,
	            grid_func & y) {
		log::debug << "Running kernel <restriction>" << std::endl;
		static_cast<child*>(this)->matvec(R, x, y);
	}


	template <class stencil>
	void residual(const stencil_op<stencil> & so,
	              const grid_func & x,
	              const grid_func & b,
	              grid_func & r) {
		log::debug << "Running kernel <residual>" << std::endl;
		static_cast<child*>(this)->residual(so, x, b, r);
	}


	void interp_add(const prolong_op & P,
	                const grid_func & coarse,
	                const grid_func & residual,
	                grid_func & fine) {
		log::debug << "Running kernel <interp_add>" << std::endl;
		static_cast<child*>(this)->interp_add(P, coarse, residual, fine);
	}


protected:
	std::shared_ptr<kernel_params> params;

};

}
#endif
