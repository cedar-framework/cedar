#include <cedar/2d/mpi/solver.h>
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

template<>
void copy_services<multilevel_wrapper<solver<nine_pt>>>(multilevel_wrapper<solver<nine_pt>> & slv, service_manager<stypes> & oldsman)
{
	service_manager<stypes> & newsman = slv.get_services();

	newsman.set_user_reg(oldsman.get_user_reg());
	newsman.set<mpi::message_passing>(oldsman.get_key<mpi::message_passing>());
	newsman.set<mpi::halo_exchange>(oldsman.get_key<mpi::halo_exchange>());
}


}}}
