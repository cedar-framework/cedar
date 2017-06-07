#include <cedar/3d/mpi/redist_solver.h>
#include <cedar/3d/cg/setup_cg_redist.h>


namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;
	template<>
	void setup_cg_redist(const kernel_params & params,
	                     const mpi::stencil_op<xxvii_pt> & so,
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

	template<>
	void setup_cg_redist(const kernel_params & params,
	                     const mpi::stencil_op<seven_pt> & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<mpi::redist_solver> * slv,
	                     std::vector<int> & nblocksv)
	{
		log::error << "Redistributed cg-solver does not support seven point stencils" << std::endl;
	}
}

}}}
