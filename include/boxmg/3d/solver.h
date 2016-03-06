#ifndef BOXMG_3D_SOLVER_H
#define BOXMG_3D_SOLVER_H

#include <array>

#include <boxmg/multilevel.h>
#include <boxmg/level.h>
#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/stencil_op.h>
#include <boxmg/3d/relax_stencil.h>
#include <boxmg/3d/inter/prolong_op.h>
#include <boxmg/3d/inter/restrict_op.h>
#include <boxmg/3d/kernel/registry.h>


namespace boxmg { namespace bmg3 {

struct bmg_level : Level<grid_func>
{
bmg_level(stencil_op&& A) :
	A(std::move(A)), P(inter::prolong_op()), SOR({{relax_stencil(), relax_stencil()}}) { R.associate(&P); }
bmg_level(stencil_op&& A, inter::prolong_op&& P) :
	A(std::move(A)), P(std::move(P)), SOR({{relax_stencil(), relax_stencil()}}) { R.associate(&P); }
	stencil_op A;
	inter::prolong_op P;
	inter::restrict_op R;
	grid_func x;
	grid_func res;
	grid_func b;
	std::array<relax_stencil, 2> SOR;
};


class solver : public multilevel<bmg_level, stencil_op, grid_func, kernel::registry>
{
public:
	solver(stencil_op&& fop);
	~solver();
	virtual int compute_num_levels(stencil_op & fop);
	virtual void setup_space(int nlevels);
};


}}

#endif

