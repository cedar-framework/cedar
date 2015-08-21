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


	void bmg2_operator_test(bmg2_operator op)
	{
		using namespace boxmg;

		log::error << "Test Error" << std::endl;
	}


	void bmg2_operator_destroy(bmg2_operator op)
	{
		using namespace boxmg::bmg2d::core;
		delete reinterpret_cast<mpi::StencilOp*>(op);
	}
}
