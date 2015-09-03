#include "solver/mpi/boxmg.h"
#include "core/mpi/stencil_op.h"
#include "boxmg-common.h"

#include "solver.h"

extern "C"
{
	bmg2_solver bmg2_solver_create(bmg2_operator op)
	{
		using namespace boxmg::bmg2d;

		core::mpi::StencilOp *sop = reinterpret_cast<core::mpi::StencilOp*>(op);

		solver::mpi::BoxMG *bmg = new solver::mpi::BoxMG(std::move(*sop));

		return reinterpret_cast<bmg2_solver>(bmg);
	}


	void bmg2_solver_run(bmg2_solver op, double **x, double *b)
	{
		using namespace boxmg::bmg2d;

		auto *bmg = reinterpret_cast<solver::mpi::BoxMG*>(op);
		auto grid = bmg->level(0).A.grid_ptr();

		core::mpi::GridFunc rhs(grid);

		auto sol = bmg->solve(rhs);
	}

	void bmg2_solver_destroy(bmg2_solver bmg)
	{
		using namespace boxmg::bmg2d::solver;

		delete reinterpret_cast<mpi::BoxMG*>(bmg);
	}

}
