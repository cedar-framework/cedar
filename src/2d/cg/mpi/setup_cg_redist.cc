#include <cedar/2d/cg/setup_cg_redist.h>


namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr2::mpi;
	void setup_cg_redist(const kernel_params & params, const mpi::stencil_op & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<mpi::redist_solver> * slv,
	                     std::vector<int> & nblocksv)
	{
		std::array<int, 2> nblocks;
		nblocks[0] = nblocksv[0];
		nblocks[1] = nblocksv[1];
		auto ret = std::make_shared<mpi::redist_solver>(so, conf, nblocks);

		*slv = ret;
	}
}

}}}
