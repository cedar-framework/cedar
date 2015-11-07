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
BoxMGLevel(StencilOp&& A, inter::ProlongOp&& P) : /*Level(A,P),*/
	A(std::move(A)), P(std::move(P)), SOR({{RelaxStencil(), RelaxStencil()}}) { R.associate(&P); }
	StencilOp    A;
	inter::ProlongOp   P;
	inter::RestrictOp  R;
	GridFunc     x;
	GridFunc     res;
	GridFunc     b;
	std::array<RelaxStencil, 2> SOR;
};

class BoxMG : public MultiLevel<BoxMGLevel>
{
public:
	BoxMG(StencilOp&& fop);
	~BoxMG() {delete[] bbd;};
	int compute_num_levels(StencilOp & fop);
	void add_level(StencilOp& fop, int num_levels);
	std::shared_ptr<kernel::Registry> kernel_registry();

private:
	GridFunc ABD;
	real_t *bbd;
};

}}}

#endif
