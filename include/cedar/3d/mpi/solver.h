#ifndef CEDAR_3D_SOLVER_MPI_CEDAR_H
#define CEDAR_3D_SOLVER_MPI_CEDAR_H

#include <memory>
#include <mpi.h>
#include <array>

#include "cedar/multilevel.h"
#include "cedar/level.h"
#include "cedar/3d/mpi/stencil_op.h"
#include "cedar/3d/relax_stencil.h"
#include "cedar/3d/inter/mpi/prolong_op.h"
#include "cedar/3d/inter/mpi/restrict_op.h"
#include "cedar/3d/kernel/mpi/registry.h"

namespace cedar { namespace cdr3 { namespace mpi {

struct BoxMGLevel : Level<cdr3::mpi::grid_func>
{
BoxMGLevel(cdr3::mpi::stencil_op&& A) : /*Level(A,P),*/
	A(std::move(A)), P(inter::mpi::prolong_op()), SOR({{relax_stencil(),relax_stencil()}}) { R.associate(&P); }
BoxMGLevel(cdr3::mpi::stencil_op&& A, inter::mpi::prolong_op&& P) : /*Level(A,P),*/
	A(std::move(A)), P(std::move(P)), SOR({{relax_stencil(),relax_stencil()}}) { R.associate(&P); }
	cdr3::mpi::grid_func     x;
	cdr3::mpi::grid_func     b;
	cdr3::mpi::grid_func     res;
	cdr3::mpi::stencil_op    A;
	inter::mpi::prolong_op   P;
	inter::mpi::restrict_op  R;
	std::array<relax_stencil,2> SOR;
};

class solver: public multilevel<BoxMGLevel,cdr3::mpi::stencil_op,cdr3::mpi::grid_func, kernel::mpi::registry>
{
public:
	solver(cdr3::mpi::stencil_op&& fop);
	solver(cdr3::mpi::stencil_op&& fop,
	       std::shared_ptr<config::reader> conf);
	~solver() {if (cg_solver_lu) bbd = new real_t[1];}
	virtual int compute_num_levels(cdr3::mpi::stencil_op & fop) override;
	MPI_Comm comm;
	virtual cdr3::mpi::grid_func solve(const cdr3::mpi::grid_func &b) override;
	virtual void solve(const cdr3::mpi::grid_func &b, cdr3::mpi::grid_func &x) override;
	virtual void setup_space(int nlevels) override;
	virtual void setup_cg_solve() override;
	void setup_halo();

private:
	bool cg_solver_lu;
	void *halo_ctx;
};

}}}

#endif
