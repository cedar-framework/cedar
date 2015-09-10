#include "operator.h"

#include "util/topo.h"
#include "core/mpi/stencil_op.h"
#include "kernel/registry.h"

#include "boxmg-common.h"

extern "C"
{
	bmg2_operator bmg2_operator_create(bmg2_topo topo)
	{
		using namespace boxmg::bmg2d;
		using namespace boxmg::bmg2d::core;

		auto & grid = *(reinterpret_cast<std::shared_ptr<mpi::GridTopo>*>(topo));

		mpi::StencilOp *sop = new mpi::StencilOp(grid);
		auto & sten = sop->stencil();
		sten.five_pt() = true;

		return reinterpret_cast<bmg2_operator>(sop);
	}


	void bmg2_operator_set(bmg2_operator op, unsigned int nvals, grid_coord coords[], double vals[])
	{
		using namespace boxmg;
		using namespace boxmg::bmg2d::core;

		auto & sten = reinterpret_cast<mpi::StencilOp*>(op)->stencil();
		auto & grid = reinterpret_cast<mpi::StencilOp*>(op)->grid();

		for (auto i: range(nvals)) {
			len_t ci = static_cast<len_t>(coords[i].i - grid.is(0) + 2);
			len_t cj = static_cast<len_t>(coords[i].j - grid.is(1) + 2);
			// boxmg likes positive stencil coefficients
			if (coords[i].dir != 0) vals[i] = -1*vals[i];
			// if (coords[i].dir == 0)
			// log::error << ci << " " << cj << " -> " << coords[i].dir << " => " << vals[i]<< std::endl;
			sten(static_cast<len_t>(ci),
			     static_cast<len_t>(cj),
			     static_cast<Dir>(coords[i].dir)) = vals[i];
		}
	}


	void bmg2_operator_apply(bmg2_operator op, const double *x, double *b)
	{
		using namespace boxmg::bmg2d;
		using namespace boxmg::bmg2d::core;

		auto *sop = reinterpret_cast<mpi::StencilOp*>(op);
		auto grid = sop->grid_ptr();

		mpi::GridFunc xgf(grid);
		int idx = 0;
		for (auto j : xgf.range(1)) {
			for (auto i : xgf.range(0)) {
				xgf(i,j) = x[idx];
				idx++;
			}
		}

		std::shared_ptr<kernel::Registry> kreg;
		kreg = std::dynamic_pointer_cast<kernel::Registry>(sop->get_registry());
		xgf.halo_ctx = sop->halo_ctx;
		kreg->halo_exchange(xgf);

		mpi::GridFunc bgf(grid);
		sop->apply(xgf, bgf);

		idx = 0;
		for (auto j : bgf.range(1)) {
			for (auto i : bgf.range(0)) {
				b[idx] = bgf(i,j);
				idx++;
			}
		}
	}


	void bmg2_operator_dump(bmg2_operator op)
	{
		using namespace boxmg::bmg2d::core;
		std::ofstream ofile;

		auto &grid = reinterpret_cast<mpi::StencilOp*>(op)->grid();

		ofile.open("op" + std::to_string(grid.coord(0)) + "-" + std::to_string(grid.coord(1)) + ".txt", std::ios::out | std::ios::trunc | std::ios::binary);
		ofile << *(reinterpret_cast<mpi::StencilOp*>(op));
		ofile.close();

	}


	void bmg2_operator_destroy(bmg2_operator op)
	{
		using namespace boxmg::bmg2d::core;

		delete reinterpret_cast<mpi::StencilOp*>(op);
	}
}
