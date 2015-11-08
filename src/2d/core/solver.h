#ifndef BOXMG_2D_CORE_SOLVER_H
#define BOXMG_2D_CORE_SOLVER_H

#include <array>

#include "multilevel.h"
#include "level.h"
#include "core/stencil_op.h"
#include "core/relax_stencil.h"
#include "inter/prolong_op.h"
#include "inter/restrict_op.h"
#include "kernel/registry.h"

namespace boxmg { namespace bmg2d {

struct BoxMGLevel : Level
{
BoxMGLevel(stencil_op&& A, inter::prolong_op&& P) : /*Level(A,P),*/
	A(std::move(A)), P(std::move(P)), SOR({{relax_stencil(), relax_stencil()}}) { R.associate(&P); }
	stencil_op    A;
	inter::prolong_op   P;
	inter::restrict_op  R;
	grid_func     x;
	grid_func     res;
	grid_func     b;
	std::array<relax_stencil, 2> SOR;
};

class solver: public MultiLevel<BoxMGLevel>
{
public:
	solver(stencil_op&& fop);
	~solver() {delete[] bbd;};
	int compute_num_levels(stencil_op & fop);
	void add_level(stencil_op& fop, int num_levels);
	std::shared_ptr<kernel::Registry> kernel_registry();

private:
	grid_func ABD;
	real_t *bbd;
};

}}

#endif
