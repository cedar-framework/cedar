#ifndef CEDAR_2D_SOLVER_MPI_CEDAR_H
#define CEDAR_2D_SOLVER_MPI_CEDAR_H

#include <memory>
#include <mpi.h>
#include <array>

#include "cedar/multilevel.h"
#include "cedar/level.h"
#include "cedar/2d/mpi/stencil_op.h"
#include "cedar/2d/relax_stencil.h"
#include "cedar/2d/inter/mpi/prolong_op.h"
#include "cedar/2d/inter/mpi/restrict_op.h"
#include "cedar/2d/kernel/mpi/registry.h"

namespace cedar { namespace cdr2 { namespace mpi {

struct boxmg_level : Level<cdr2::mpi::grid_func>
{
boxmg_level(cdr2::mpi::stencil_op&& A) : /*Level(A,P),*/
	A(std::move(A)), P(inter::mpi::prolong_op()), SOR({{relax_stencil(),relax_stencil()}}) { R.associate(&P); }
boxmg_level(cdr2::mpi::stencil_op&& A, inter::mpi::prolong_op&& P) : /*Level(A,P),*/
	A(std::move(A)), P(std::move(P)), SOR({{relax_stencil(),relax_stencil()}}) { R.associate(&P); }
	cdr2::mpi::grid_func     x;
	cdr2::mpi::grid_func     b;
	cdr2::mpi::grid_func     res;
	cdr2::mpi::stencil_op    A;
	inter::mpi::prolong_op   P;
	inter::mpi::restrict_op  R;
	std::array<relax_stencil,2> SOR;
};

class solver: public multilevel<boxmg_level,cdr2::mpi::stencil_op, cdr2::mpi::grid_func, kernel::mpi::registry>
{
public:
	solver(cdr2::mpi::stencil_op&& fop);
	solver(cdr2::mpi::stencil_op&& fop,
	       std::shared_ptr<config::reader> conf);
	~solver() {if (cg_solver_lu) bbd = new real_t[1];}
	virtual int compute_num_levels(cdr2::mpi::stencil_op & fop) override;
	MPI_Comm comm;
	virtual cdr2::mpi::grid_func solve(const cdr2::mpi::grid_func &b) override;
	virtual void solve(const cdr2::mpi::grid_func &b, cdr2::mpi::grid_func &x) override;
	virtual void setup_cg_solve() override;
	virtual void setup_space(int nlevels) override;
	void setup_halo();

private:
	bool cg_solver_lu;
	void *halo_ctx;
};

}}}

#endif
