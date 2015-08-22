#include "operator.h"

#include "util/topo.h"
#include "core/mpi/stencil_op.h"

#include "boxmg-common.h"

extern "C"
{
	bmg2_operator bmg2_operator_create(unsigned int nx, unsigned int ny)
	{
		using namespace boxmg::bmg2d;
		using namespace boxmg::bmg2d::core;
		auto grid = util::create_topo(MPI_COMM_WORLD, nx, ny);
		return reinterpret_cast<bmg2_operator>(new mpi::StencilOp(grid));
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
		}
	}


	void bmg2_operator_destroy(bmg2_operator op)
	{
		using namespace boxmg::bmg2d::core;
		delete reinterpret_cast<mpi::StencilOp*>(op);
	}
}
