#include <boxmg/2d/cg/setup_cg_redist.h>


namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg2d::mpi;
	void setup_cg_redist(const mpi::stencil_op & so,
	                     std::shared_ptr<mpi::solver> * slv)
	{
	}
}

}}}
