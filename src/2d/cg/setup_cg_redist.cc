#include <boxmg/2d/cg/setup_cg_redist.h>


namespace boxmg { namespace bmg2d { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg2d::mpi;
	void setup_cg_redist(const mpi::stencil_op & so,
	                     std::shared_ptr<mpi::redist_solver> * slv)
	{
		log::status << "Setting up cg redist" << std::endl;
		std::array<int, 2> nblocks({2,2});
		auto ret = std::make_shared<mpi::redist_solver>(so, nblocks);

		*slv = ret;
	}
}

}}}
