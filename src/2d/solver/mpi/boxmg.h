#ifndef BOXMG_2D_SOLVER_MPI_BOXMG_H
#define BOXMG_2D_SOLVER_MPI_BOXMG_H

#include <memory>
#include <mpi.h>
#include <array>

#include "../multilevel.h"
#include "core/mpi/stencil_op.h"
#include "core/relax_stencil.h"
#include "inter/mpi/prolong_op.h"
#include "inter/mpi/restrict_op.h"
#include "kernel/registry.h"

namespace boxmg { namespace bmg2d { namespace solver { namespace mpi {

struct BoxMGLevel : Level
{
BoxMGLevel(bmg2d::mpi::StencilOp&& A, inter::mpi::ProlongOp&& P) : /*Level(A,P),*/
	A(std::move(A)), P(std::move(P)), SOR({{RelaxStencil(),RelaxStencil()}}) { R.associate(&P); }
	bmg2d::mpi::GridFunc     x;
	bmg2d::mpi::GridFunc     b;
	bmg2d::mpi::GridFunc     res;
	bmg2d::mpi::StencilOp    A;
	inter::mpi::ProlongOp   P;
	inter::mpi::RestrictOp  R;
	std::array<RelaxStencil,2> SOR;
};

class BoxMG : public MultiLevel<BoxMGLevel,bmg2d::mpi::GridFunc>
{
public:
	BoxMG(bmg2d::mpi::StencilOp&& fop);
	~BoxMG() {if (cg_solver_lu) delete[] bbd;};
	int compute_num_levels(bmg2d::mpi::StencilOp & fop);
	void add_level(bmg2d::mpi::StencilOp& fop, int num_levels);
	MPI_Comm comm;
	std::shared_ptr<kernel::Registry> kernel_registry();
	virtual bmg2d::mpi::GridFunc solve(const bmg2d::mpi::GridFunc &b);
	virtual void solve(const bmg2d::mpi::GridFunc &b, bmg2d::mpi::GridFunc &x);

private:
	GridFunc ABD;
	bool cg_solver_lu;
	real_t *bbd;
	void *halo_ctx;
};

}}}}

#endif
