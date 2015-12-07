#ifndef BOXMG_2D_KERNEL_REGISTRY_H
#define BOXMG_2D_KERNEL_REGISTRY_H


#include "boxmg/types.h"
#include "boxmg/cycle/types.h"

#include "boxmg/2d/stencil_op.h"
#include "boxmg/2d/inter/prolong_op.h"
#include "boxmg/2d/relax_stencil.h"
#include "boxmg/2d/grid_func.h"
#include "boxmg/2d/mpi/grid_topo.h"
#include "boxmg/2d/mpi/stencil_op.h"
#include "boxmg/2d/mpi/grid_func.h"


namespace boxmg { namespace bmg2d {
			class solver;
}}


namespace boxmg { namespace bmg2d { namespace kernel {

class Registry : public KernelRegistry
{
public:
	void setup_interp(int kf, int kc, int nog, const stencil_op & fop,
	                  const stencil_op &cop, inter::prolong_op & P);

	void galerkin_prod(int kf, int kc, int nog,
	                   const inter::prolong_op & P,
	                   const stencil_op & fop,
	                   stencil_op & cop);

	void setup_relax(const stencil_op & so,
	                 relax_stencil & sor);

	void setup_relax_x(const stencil_op & so,
	                   relax_stencil & sor);

	void setup_relax_y(const stencil_op & so,
	                   relax_stencil & sor);

	void setup_cg_lu(const stencil_op & so,
	                 grid_func & ABD);

	void relax(const stencil_op & so,
	           grid_func & x,
	           const grid_func & b,
	           const relax_stencil & sor,
	           cycle::Dir cycle_dir);

	void relax_lines_x(const stencil_op & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func &res,
	                   cycle::Dir cycle_dir);

	void relax_lines_y(const stencil_op & so,
	                   grid_func & x,
	                   const grid_func & b,
	                   const relax_stencil & sor,
	                   grid_func &res,
	                   cycle::Dir cycle_dir);

	void solve_cg(grid_func &x,
	              const grid_func &b,
	              const grid_func &ABD,
	              real_t *bbd);

	void setup_nog(mpi::grid_topo &topo,
	               len_t min_coarse, int *nog);

	void halo_setup(mpi::grid_topo &topo,
	                void **halo_ctx);
	void halo_exchange(mpi::grid_func &f);
	void halo_exchange(const mpi::grid_func &f, void *halo_ctx);
	void halo_stencil_exchange(mpi::stencil_op & so);
	void setup_cg_boxmg(const stencil_op & so,
	                    std::shared_ptr<solver> *bmg);
	void solve_cg_boxmg(const solver &bmg,
	                    grid_func &x,
	                    const grid_func &b);
	void matvec(const stencil_op & so,
	            const grid_func &x,
	            grid_func &b);
};

}}}


#endif
