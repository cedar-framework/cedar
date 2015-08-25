#include "operator.h"

#include "util/topo.h"
#include "core/mpi/stencil_op.h"

#include "boxmg-common.h"

extern "C"
{
	bmg2_operator bmg2_operator_create(bmg2_topo topo)
	{
		using namespace boxmg::bmg2d;
		using namespace boxmg::bmg2d::core;

		auto grid = *(reinterpret_cast<std::shared_ptr<mpi::GridTopo>*>(topo));

		mpi::StencilOp *sop = new mpi::StencilOp(grid);
		auto & sten = sop->stencil();
		sten.five_pt() = true;
		return reinterpret_cast<bmg2_operator>(sop);
	}


	void bmg2_operator_set(bmg2_operator op, unsigned int nvals, grid_coord coords[], double vals[])
	{
		using namespace boxmg;
		using namespace boxmg::bmg2d::core;

		GridStencil & sten = reinterpret_cast<mpi::StencilOp*>(op)->stencil();

		for (auto i: range(nvals)) {
			sten(static_cast<len_t>(coords[i].i),
			     static_cast<len_t>(coords[i].j),
			     static_cast<Dir>(coords[i].dir)) = vals[i];
			log::error << coords[i].i << " " << coords[i].j << " -> " << coords[i].dir << std::endl;
		}
	}


	void bmg2_operator_destroy(bmg2_operator op)
	{
		using namespace boxmg::bmg2d::core;

		std::ofstream ofile;

		ofile.open("op.txt", std::ios::out | std::ios::trunc | std::ios::binary);
		ofile << *(reinterpret_cast<mpi::StencilOp*>(op));
		ofile.close();

		delete reinterpret_cast<mpi::StencilOp*>(op);
	}
}
