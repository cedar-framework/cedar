#include <boxmg/types.h>
#include <boxmg/util/log.h>
#include <boxmg/3d/mpi/solver.h>
#include <boxmg/3d/mpi/stencil_op.h>

#include <boxmg/3d/interface/c/solver.h>

extern "C"
{
	bmg3_solver bmg3_solver_create(bmg3_operator *op)
	{
		using namespace boxmg::bmg3;

		mpi::stencil_op *sop = reinterpret_cast<mpi::stencil_op*>(*op);

		boxmg::log::info << "Beginning setup phase" << std::endl;
		mpi::solver *bmg = new mpi::solver(std::move(*sop));
		boxmg::log::info << "Setup phase complete" << std::endl;
		bmg->level(-1).x = mpi::grid_func::zeros_like(bmg->level(-1).res);
		bmg->level(-1).b = mpi::grid_func::zeros_like(bmg->level(-1).res);
		*op = reinterpret_cast<bmg3_operator>(&bmg->level(-1).A);

		return reinterpret_cast<bmg3_solver>(bmg);
	}


	void bmg3_solver_run(bmg3_solver op, double *x, const double *b)
	{
		using namespace boxmg::bmg3;

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

	void bmg3_solver_destroy(bmg3_solver bmg)
	{
		using namespace boxmg::bmg3;

		delete reinterpret_cast<mpi::solver*>(bmg);
	}

}
