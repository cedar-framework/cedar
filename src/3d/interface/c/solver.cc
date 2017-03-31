#include <cedar/types.h>
#include <cedar/util/log.h>
#include <cedar/3d/mpi/solver.h>
#include <cedar/3d/mpi/stencil_op.h>

#include <cedar/3d/interface/c/solver.h>

extern "C"
{
	cdr3_solver cdr3_solver_create(cdr3_operator *op)
	{
		using namespace cedar::cdr3;

		mpi::stencil_op *sop = reinterpret_cast<mpi::stencil_op*>(*op);

		cedar::log::info << "Beginning setup phase" << std::endl;
		mpi::solver *bmg = new mpi::solver(std::move(*sop));
		cedar::log::info << "Setup phase complete" << std::endl;
		bmg->level(-1).x = mpi::grid_func::zeros_like(bmg->level(-1).res);
		bmg->level(-1).b = mpi::grid_func::zeros_like(bmg->level(-1).res);
		*op = reinterpret_cast<cdr3_operator>(&bmg->level(-1).A);

		return reinterpret_cast<cdr3_solver>(bmg);
	}


	void cdr3_solver_run(cdr3_solver op, double *x, const double *b)
	{
		using namespace cedar::cdr3;

		auto *bmg = reinterpret_cast<mpi::solver*>(op);
		auto grid = bmg->level(-1).A.grid_ptr();

		auto & rhs = bmg->level(-1).b;

		int idx = 0;
		for (auto k : rhs.range(2)) {
			for (auto j : rhs.range(1)) {
				for (auto i : rhs.range(0)) {
					rhs(i,j,k) = b[idx];
					idx++;
				}
			}
		}

		auto & sol = bmg->level(-1).x;
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

	void cdr3_solver_destroy(cdr3_solver bmg)
	{
		using namespace cedar::cdr3;

		delete reinterpret_cast<mpi::solver*>(bmg);
	}

}
