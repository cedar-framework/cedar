#ifndef BOXMG_2D_SOLVER_MPI_BOXMG_H
#define BOXMG_2D_SOLVER_MPI_BOXMG_H

#include <memory>
#include <mpi.h>
#include <array>

#include "boxmg/multilevel.h"
#include "boxmg/level.h"
#include "boxmg/2d/mpi/stencil_op.h"
#include "boxmg/2d/relax_stencil.h"
#include "boxmg/2d/inter/mpi/prolong_op.h"
#include "boxmg/2d/inter/mpi/restrict_op.h"
#include "boxmg/2d/kernel/mpi/registry.h"

namespace boxmg { namespace bmg2d { namespace mpi {

struct boxmg_level : Level<bmg2d::mpi::grid_func>
{
boxmg_level(bmg2d::mpi::stencil_op&& A) : /*Level(A,P),*/
	A(std::move(A)), P(inter::mpi::prolong_op()), SOR({{relax_stencil(),relax_stencil()}}) { R.associate(&P); }
boxmg_level(bmg2d::mpi::stencil_op&& A, inter::mpi::prolong_op&& P) : /*Level(A,P),*/
	A(std::move(A)), P(std::move(P)), SOR({{relax_stencil(),relax_stencil()}}) { R.associate(&P); }
	bmg2d::mpi::grid_func     x;
	bmg2d::mpi::grid_func     b;
	bmg2d::mpi::grid_func     res;
	bmg2d::mpi::stencil_op    A;
	inter::mpi::prolong_op   P;
	inter::mpi::restrict_op  R;
	std::array<relax_stencil,2> SOR;
};

class solver: public multilevel<boxmg_level,bmg2d::mpi::stencil_op, bmg2d::mpi::grid_func, kernel::mpi::registry>
{
public:
	solver(bmg2d::mpi::stencil_op&& fop);
	~solver() {if (cg_solver_lu) bbd = new real_t[1];}
	virtual int compute_num_levels(bmg2d::mpi::stencil_op & fop);
	MPI_Comm comm;
	virtual bmg2d::mpi::grid_func solve(const bmg2d::mpi::grid_func &b);
	virtual void solve(const bmg2d::mpi::grid_func &b, bmg2d::mpi::grid_func &x);
	virtual void setup_cg_solve();
	virtual void setup_space(int nlevels);
	void setup_halo();

private:
	bool cg_solver_lu;
	void *halo_ctx;
};

}}}

#endif
