#ifndef BOXMG_2D_SOLVER_BOXMG_H
#define BOXMG_2D_SOLVER_BOXMG_H

#include <array>

#include "multilevel.h"
#include "core/stencil_op.h"
#include "core/relax_stencil.h"
#include "inter/prolong_op.h"
#include "inter/restrict_op.h"
#include "kernel/registry.h"

namespace boxmg { namespace bmg2d { namespace solver {

struct BoxMGLevel : Level
{
BoxMGLevel(core::StencilOp&& A, inter::ProlongOp&& P) : /*Level(A,P),*/
	A(std::move(A)), P(std::move(P)), SOR({{core::RelaxStencil(), core::RelaxStencil()}}) { R.associate(&P); }
	core::StencilOp    A;
	inter::ProlongOp   P;
	inter::RestrictOp  R;
	std::array<core::RelaxStencil, 2> SOR;
};

class BoxMG : public MultiLevel<BoxMGLevel>
{
public:
	BoxMG(core::StencilOp&& fop);
	~BoxMG() {};
	int compute_num_levels(core::StencilOp & fop);
	void add_level(core::StencilOp& fop, int num_levels);
	std::shared_ptr<kernel::Registry> kernel_registry();

private:
	core::GridFunc ABD;
};

}}}

#endif
