#include <boxmg/3d/cg/setup_cg_redist.h>


namespace boxmg { namespace bmg3 { namespace kernel {

namespace impls
{
	namespace mpi = boxmg::bmg3::mpi;
	void setup_cg_redist(const mpi::stencil_op & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<mpi::redist_solver> * slv,
	                     std::vector<int> & nblocksv)
	{
		std::array<int, 3> nblocks;
		for (auto i : range(3)) {
			nblocks[i] = nblocksv[i];
		}

		auto ret = std::make_shared<mpi::redist_solver>(so, conf, nblocks);

		*slv = ret;
	}
}

}}}
