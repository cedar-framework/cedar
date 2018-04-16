#include <cedar/2d/mpi/redist_solver.h>

namespace cedar { namespace cdr2 { namespace mpi {

template<>
	std::unique_ptr<cdr2::stencil_op<nine_pt>> create_operator(topo_ptr topo)
{
	return std::make_unique<cdr2::stencil_op<nine_pt>>(topo->nlocal(0)-2, topo->nlocal(1)-2);
}

template<>
	std::unique_ptr<mpi::stencil_op<nine_pt>> create_operator(topo_ptr topo)
{
	return std::make_unique<mpi::stencil_op<nine_pt>>(topo);
}


}}}
