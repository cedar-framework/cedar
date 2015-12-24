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

struct BoxMGLevel : Level<bmg2d::mpi::grid_func>
{
BoxMGLevel(bmg2d::mpi::stencil_op&& A, inter::mpi::prolong_op&& P) : /*Level(A,P),*/
	A(std::move(A)), P(std::move(P)), SOR({{relax_stencil(),relax_stencil()}}) { R.associate(&P); }
	bmg2d::mpi::grid_func     x;
	bmg2d::mpi::grid_func     b;
	bmg2d::mpi::grid_func     res;
	bmg2d::mpi::stencil_op    A;
	inter::mpi::prolong_op   P;
	inter::mpi::restrict_op  R;
	std::array<relax_stencil,2> SOR;
};

class solver: public multilevel<BoxMGLevel,bmg2d::mpi::grid_func, kernel::mpi::registry>
{
public:
	solver(bmg2d::mpi::stencil_op&& fop);
	~solver() {if (cg_solver_lu) delete[] bbd;};
	int compute_num_levels(bmg2d::mpi::stencil_op & fop);
	void add_level(bmg2d::mpi::stencil_op& fop, int num_levels);
	MPI_Comm comm;
	std::shared_ptr<kernel::mpi::registry> kernel_registry();
	virtual bmg2d::mpi::grid_func solve(const bmg2d::mpi::grid_func &b);
	virtual void solve(const bmg2d::mpi::grid_func &b, bmg2d::mpi::grid_func &x);

private:
	bmg2d::mpi::grid_func ABD;
	bool cg_solver_lu;
	real_t *bbd;
	void *halo_ctx;
};

}}}

#endif
