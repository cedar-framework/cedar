#ifndef BOXMG_KERNEL_REGISTRY_H
#define BOXMG_KERNEL_REGISTRY_H

#include <memory>
#include <map>
#include "kernel_manager.h"
#include <boxmg/cycle/types.h>
#include <boxmg/kernel_name.h>
#include <boxmg/types.h>

namespace boxmg {


template <class stencil_op, class relax_stencil, class prolong_op, class grid_func>

class kernel_registry
{
public:
	template <typename S, typename T>
		void add(const std::string & kname, S&& kid, const T&cmd)
	{
		avail[kname].add(std::forward<S>(kid), cmd);
	}


	virtual void set(const std::string & kname, const std::string & kid)
	{
		active[kname] = avail[kname].at(kid);
	}


	template <typename... Args>
		void run(const std::string & kname, Args&&... args)
	{
		active.run(kname, std::forward<decltype(args)>(args)...);
	}

	void setup_interp(int kf, int kc, int nog, const stencil_op & fop,
	                  const stencil_op &cop, prolong_op & P)
	{
		active.run(kernel_name::setup_interp,
		           static_cast<int>(kf),
		           static_cast<int>(kc),
		           static_cast<int>(nog),
		           fop,cop,P);
	}


	void galerkin_prod(int kf, int kc, int nog,
	                   const prolong_op & P,
	                   const stencil_op & fop,
	                   stencil_op & cop)
	{
		active.run(kernel_name::galerkin_prod,
		           static_cast<int>(kf),
		           static_cast<int>(kc),
		           static_cast<int>(nog),
		           P,fop,cop);
	}


	void setup_relax(const stencil_op & so,
	                 relax_stencil & sor)
	{
		active.run(kernel_name::setup_relax, so, sor);
	}


	void setup_relax_x(const stencil_op & so,
	                   relax_stencil & sor)
	{
		active.run(kernel_name::setup_relax_x, so, sor);
	}


	void setup_relax_y(const stencil_op & so,
	                   relax_stencil & sor)
	{
		active.run(kernel_name::setup_relax_y, so, sor);
	}


	void setup_cg_lu(const stencil_op & so,
	                 grid_func & ABD)
	{
		active.run(kernel_name::setup_cg_lu, so, ABD);
	}


	void relax(const stencil_op & so,
	           grid_func & x,
	           const grid_func & b,
	           const relax_stencil & sor,
	           cycle::Dir cdir)
	{
		active.run(kernel_name::relax, so, x, b, sor, static_cast<cycle::Dir>(cdir));
	}


	void relax_lines_x(const stencil_op & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func &res,
	                   cycle::Dir cdir)
	{
		active.run(kernel_name::relax_lines_x, so, x, b, sor, res, static_cast<cycle::Dir>(cdir));
	}


	void relax_lines_y(const stencil_op & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func &res,
	                   cycle::Dir cdir)
	{
		active.run(kernel_name::relax_lines_y, so, x, b, sor, res, static_cast<cycle::Dir>(cdir));
	}


	void solve_cg(grid_func &x,
	              const grid_func &b,
	              const grid_func &ABD,
	              real_t *bbd)
	{
		active.run(kernel_name::solve_cg, x, b, ABD, static_cast<real_t*>(bbd));
	}

protected:
	kernel_manager active;
	std::map<std::string, kernel_manager> avail;
};

}
#endif
