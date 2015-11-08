#ifndef BOXMG_2D_KERNEL_REGISTRY_H
#define BOXMG_2D_KERNEL_REGISTRY_H


#include "boxmg-common.h"

#include "core/stencil_op.h"
#include "inter/prolong_op.h"
#include "core/relax_stencil.h"
#include "core/grid_func.h"
#include "core/mpi/grid_topo.h"
#include "core/mpi/stencil_op.h"
#include "core/mpi/grid_func.h"


namespace boxmg { namespace bmg2d {
			class solver;
}}


namespace boxmg { namespace bmg2d { namespace kernel {

class Registry : public KernelRegistry
{
public:
	void setup_interp(int kf, int kc, int nog, const StencilOp & fop,
	                  const StencilOp &cop, inter::prolong_op & P);

	void galerkin_prod(int kf, int kc, int nog,
	                   const inter::prolong_op & P,
	                   const StencilOp & fop,
	                   StencilOp & cop);

	void setup_relax(const StencilOp & so,
	                 RelaxStencil & sor);

	void setup_relax_x(const StencilOp & so,
	                   RelaxStencil & sor);

	void setup_relax_y(const StencilOp & so,
	                   RelaxStencil & sor);

	void setup_cg_lu(const StencilOp & so,
	                 grid_func & ABD);

	void relax(const StencilOp & so,
	           grid_func & x,
	           const grid_func & b,
	           const RelaxStencil & sor,
	           cycle::Dir cycle_dir);

	void relax_lines_x(const StencilOp & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const RelaxStencil & sor,
	                   grid_func &res,
	                   cycle::Dir cycle_dir);

	void relax_lines_y(const StencilOp & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const RelaxStencil & sor,
	                   grid_func &res,
	                   cycle::Dir cycle_dir);

	void solve_cg(grid_func &x,
	              const grid_func &b,
	              const grid_func &ABD,
	              real_t *bbd);

	void setup_nog(mpi::GridTopo &topo,
	               len_t min_coarse, int *nog);

	void halo_setup(mpi::GridTopo &topo,
	                void **halo_ctx);
	void halo_exchange(mpi::grid_func &f);
	void halo_exchange(const mpi::grid_func &f, void *halo_ctx);
	void halo_stencil_exchange(mpi::StencilOp & so);
	void setup_cg_boxmg(const StencilOp & so,
	                    std::shared_ptr<solver> *bmg);
	void solve_cg_boxmg(const solver &bmg,
	                    grid_func &x,
	                    const grid_func &b);
	void matvec(const StencilOp & so,
	            const grid_func &x,
	            grid_func &b);
};

}}}


#endif
