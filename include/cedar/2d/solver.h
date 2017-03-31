#ifndef CEDAR_2D_CORE_SOLVER_H
#define CEDAR_2D_CORE_SOLVER_H

#include <array>

#include "cedar/multilevel.h"
#include "cedar/level.h"
#include "cedar/2d/stencil_op.h"
#include "cedar/2d/relax_stencil.h"
#include "cedar/2d/inter/prolong_op.h"
#include "cedar/2d/inter/restrict_op.h"
#include "cedar/2d/kernel/registry.h"

namespace cedar { namespace cdr2 {

struct BoxMGLevel : Level<grid_func>
{
BoxMGLevel(stencil_op&& A) : /*Level(A,P),*/
	A(std::move(A)), P(inter::prolong_op()), SOR({{relax_stencil(), relax_stencil()}}) { R.associate(&P); }
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

class solver: public multilevel<BoxMGLevel, stencil_op, grid_func, kernel::registry>
{
public:
	solver(stencil_op&& fop);
	solver(stencil_op&& fop,
	       std::shared_ptr<config::reader> conf);
	~solver();
	virtual int compute_num_levels(stencil_op & fop);
	virtual void setup_space(int nlevels);
};

}}

#endif
