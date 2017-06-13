#include <cedar/types.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/mpi/stencil_op.h>

#include <cedar/2d/interface/c/solver.h>
#include <cedar/2d/interface/c/types.h>

extern "C"
{
	bmg2_solver bmg2_solver_create(bmg2_operator *op)
	{
		using namespace cedar::cdr2;

		auto & sop = reinterpret_cast<op_container*>(*op)->op;

		auto *bmg = new mpi::solver<nine_pt>(sop);
		auto & flevel = bmg->levels.get(0);
		flevel.x = mpi::grid_func::zeros_like(flevel.res);
		flevel.b = mpi::grid_func::zeros_like(flevel.res);

		return reinterpret_cast<bmg2_solver>(bmg);
	}


	void bmg2_solver_run(bmg2_solver op, double *x, const double *b)
	{
		using namespace cedar::cdr2;

		auto *bmg = reinterpret_cast<mpi::solver<nine_pt>*>(op);

		auto & rhs = bmg->levels.get(0).b;

		int idx = 0;
		for (auto j : rhs.range(1)) {
			for (auto i : rhs.range(0)) {
				rhs(i,j) = b[idx];
				idx++;
			}
		}

		auto & sol = bmg->levels.get(0).x;
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
		using namespace cedar::cdr2;

		delete reinterpret_cast<mpi::solver<nine_pt>*>(bmg);
	}

}
