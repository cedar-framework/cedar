#include "solver/mpi/boxmg.h"
#include "core/mpi/stencil_op.h"
#include "boxmg-common.h"

#include "solver.h"

extern "C"
{
	bmg2_solver bmg2_solver_create(bmg2_operator *op)
	{
		using namespace boxmg::bmg2d;

		mpi::StencilOp *sop = reinterpret_cast<mpi::StencilOp*>(*op);

		solver::mpi::BoxMG *bmg = new solver::mpi::BoxMG(std::move(*sop));
		bmg->level(-1).x = mpi::GridFunc::zeros_like(bmg->level(-1).res);
		bmg->level(-1).b = mpi::GridFunc::zeros_like(bmg->level(-1).res);
		*op = reinterpret_cast<bmg2_operator>(&bmg->level(-1).A);

		return reinterpret_cast<bmg2_solver>(bmg);
	}


	void bmg2_solver_run(bmg2_solver op, double *x, const double *b)
	{
		using namespace boxmg::bmg2d;

		auto *bmg = reinterpret_cast<solver::mpi::BoxMG*>(op);
		auto grid = bmg->level(-1).A.grid_ptr();

		auto & rhs = bmg->level(-1).b;

		int idx = 0;
		for (auto j : rhs.range(1)) {
			for (auto i : rhs.range(0)) {
				rhs(i,j) = b[idx];
				idx++;
			}
		}

		auto & sol = bmg->level(-1).x;
		sol.set(0.0);
		bmg->solve(rhs, sol);

		idx = 0;
		for (auto j : sol.range(1)) {
			for (auto i : sol.range(0)) {
				x[idx] = sol(i,j);
				idx++;
			}
		}
	}

	void bmg2_solver_destroy(bmg2_solver bmg)
	{
		using namespace boxmg::bmg2d::solver;

		delete reinterpret_cast<mpi::BoxMG*>(bmg);
	}

}
