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


namespace boxmg { namespace bmg2d { namespace solver {
			class BoxMG;
}}}


namespace boxmg { namespace bmg2d { namespace kernel {

class Registry : public KernelRegistry
{
public:
	void setup_interp(int kf, int kc, int nog, const core::StencilOp & fop,
	                  const core::StencilOp &cop, inter::ProlongOp & P);

	void galerkin_prod(int kf, int kc, int nog,
	                   const inter::ProlongOp & P,
	                   const core::StencilOp & fop,
	                   core::StencilOp & cop);

	void setup_relax(const core::StencilOp & so,
	                 core::RelaxStencil & sor);

	void setup_cg_lu(const core::StencilOp & so,
	                 core::GridFunc & ABD);

	void relax(const core::StencilOp & so,
	           core::GridFunc & x,
	           const core::GridFunc & b,
	           const core::RelaxStencil & sor,
	           cycle::Dir cycle_dir);

	void solve_cg(core::GridFunc &x,
	              const core::GridFunc &b,
	              const core::GridFunc &ABD);

	void setup_nog(core::mpi::GridTopo &topo,
	               len_t min_coarse, int *nog);

	void halo_setup(core::mpi::GridTopo &topo,
	                void **halo_ctx);
	void halo_exchange(core::mpi::GridFunc &f);
	void halo_exchange(const core::mpi::GridFunc &f, void *halo_ctx);
	void halo_stencil_exchange(core::mpi::StencilOp & so);
	void setup_cg_boxmg(const core::StencilOp & so,
	                    std::shared_ptr<solver::BoxMG> *solver);
	void solve_cg_boxmg(const solver::BoxMG &bmg,
	                    core::GridFunc &x,
	                    const core::GridFunc &b);
};

}}}


#endif
