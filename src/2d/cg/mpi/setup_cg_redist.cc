#include <cedar/2d/mpi/redist_solver.h>
#include <cedar/2d/cg/setup_cg_redist.h>


namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr2::mpi;

	template<>
	void setup_cg_redist(const kernel_params & params,
	                     halo_exchanger *halof,
	                     const mpi::stencil_op<nine_pt> & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<mpi::redist_solver> * slv,
	                     std::vector<int> & nblocksv)
	{
		std::array<int, 2> nblocks;
		nblocks[0] = nblocksv[0];
		nblocks[1] = nblocksv[1];
		auto ret = std::make_shared<mpi::redist_solver>(so, halof, conf, nblocks);

		*slv = ret;
	}


	template<>
	void setup_cg_redist(const kernel_params & params,
	                     halo_exchanger *halof,
	                     const mpi::stencil_op<five_pt> & so,
	                     std::shared_ptr<config::reader> conf,
	                     std::shared_ptr<mpi::redist_solver> * slv,
	                     std::vector<int> & nblocksv)
	{
		log::error << "Redist cg solver is not implemented for five point stencils" << std::endl;
	}
}

}}}
