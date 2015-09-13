#ifndef BOXMG_2D_SOLVER_MPI_BOXMG_H
#define BOXMG_2D_SOLVER_MPI_BOXMG_H

#include <memory>
#include <mpi.h>
#include <array>

#include "../multilevel.h"
#include "core/mpi/stencil_op.h"
#include "core/relax_stencil.h"
#include "inter/mpi/prolong_op.h"
#include "inter/mpi/restrict_op.h"
#include "kernel/registry.h"

namespace boxmg { namespace bmg2d { namespace solver { namespace mpi {

struct BoxMGLevel : Level
{
BoxMGLevel(core::mpi::StencilOp&& A, inter::mpi::ProlongOp&& P) : /*Level(A,P),*/
	A(std::move(A)), P(std::move(P)), SOR({{core::RelaxStencil(),core::RelaxStencil()}}) { R.associate(&P); }
	core::mpi::StencilOp    A;
	inter::mpi::ProlongOp   P;
	inter::mpi::RestrictOp  R;
	std::array<core::RelaxStencil,2> SOR;
};

class BoxMG : public MultiLevel<BoxMGLevel,core::mpi::GridFunc>
{
public:
	BoxMG(core::mpi::StencilOp&& fop);
	~BoxMG() {};
	int compute_num_levels(core::mpi::StencilOp & fop);
	void add_level(core::mpi::StencilOp& fop, int num_levels);
	MPI_Comm comm;
	std::shared_ptr<kernel::Registry> kernel_registry();
	virtual core::mpi::GridFunc solve(const core::mpi::GridFunc &b);

private:
	core::GridFunc ABD;
	void *halo_ctx;
};

}}}}

#endif
