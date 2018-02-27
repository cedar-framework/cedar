#include <cedar/types.h>
#include <cedar/util/log.h>
#include <cedar/3d/mpi/solver.h>
#include <cedar/3d/mpi/stencil_op.h>

#include <cedar/3d/interface/c/solver.h>
#include <cedar/3d/interface/c/types.h>

extern "C"
{
	bmg3_solver bmg3_solver_create(bmg3_operator *op)
	{
		using namespace cedar::cdr3;

		auto & sop = reinterpret_cast<op_container*>(*op)->op;

		cedar::log::info << "Beginning setup phase" << std::endl;

		auto *bmg = new mpi::solver<xxvii_pt>(sop);
		cedar::log::info << "Setup phase complete" << std::endl;
		auto & flevel = bmg->levels.get(0);
		flevel.x = mpi::grid_func::zeros_like(flevel.res);
		flevel.b = mpi::grid_func::zeros_like(flevel.res);

		return reinterpret_cast<bmg3_solver>(bmg);
	}


	void bmg3_solver_run(bmg3_solver op, double *x, const double *b)
	{
		using namespace cedar::cdr3;

		auto *bmg = reinterpret_cast<mpi::solver<xxvii_pt>*>(op);

		auto & rhs = bmg->levels.get(0).b;

		int idx = 0;
		for (auto k : rhs.range(2)) {
			for (auto j : rhs.range(1)) {
				for (auto i : rhs.range(0)) {
					rhs(i,j,k) = b[idx];
					idx++;
				}
			}
		}

		auto & sol = bmg->levels.get(0).x;
		sol.set(0.0);
		bmg->solve(rhs, sol);

		idx = 0;
		for (auto k : sol.range(2)) {
			for (auto j : sol.range(1)) {
				for (auto i : sol.range(0)) {
					x[idx] = sol(i,j,k);
					idx++;
				}
			}
		}
	}

	void bmg3_solver_destroy(bmg3_solver bmg)
	{
		using namespace cedar::cdr3;

		delete reinterpret_cast<mpi::solver<xxvii_pt>*>(bmg);
	}

}
