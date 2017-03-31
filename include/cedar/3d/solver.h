#ifndef CEDAR_3D_SOLVER_H
#define CEDAR_3D_SOLVER_H

#include <array>

#include <cedar/multilevel.h>
#include <cedar/level.h>
#include <cedar/2d/solver.h>
#include <cedar/3d/grid_func.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/relax_stencil.h>
#include <cedar/3d/inter/prolong_op.h>
#include <cedar/3d/inter/restrict_op.h>
#include <cedar/3d/kernel/registry.h>


namespace cedar { namespace cdr3 {

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
	std::vector<std::unique_ptr<::cedar::cdr2::solver>> planes_xy;
	std::vector<std::unique_ptr<::cedar::cdr2::solver>> planes_xz;
	std::vector<std::unique_ptr<::cedar::cdr2::solver>> planes_yz;
};


class solver : public multilevel<bmg_level, stencil_op, grid_func, kernel::registry>
{
public:
	solver(stencil_op&& fop);
	solver(stencil_op&& fop,
	       std::shared_ptr<config::reader> conf);
	~solver();
	virtual int compute_num_levels(stencil_op & fop) override;
	virtual void setup_relax_plane(stencil_op & sop, bmg_level & level) override;
	virtual void relax_plane(const stencil_op & so, grid_func & x,
	                         const grid_func & b, cycle::Dir cdir,
	                         bmg_level & level) override;
	virtual void setup_space(int nlevels) override;
};


}}

#endif

