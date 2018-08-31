#include <cedar/3d/mpi/solver.h>
#include <cedar/3d/mpi/redist_solver.h>


namespace cedar { namespace cdr3 { namespace mpi {

template<>
	std::unique_ptr<cdr3::stencil_op<xxvii_pt>> create_operator(topo_ptr topo)
{
	return std::make_unique<cdr3::stencil_op<xxvii_pt>>(topo->nlocal(0)-2, topo->nlocal(1)-2, topo->nlocal(2)-2);
}

template<>
	std::unique_ptr<mpi::stencil_op<xxvii_pt>> create_operator(topo_ptr topo)
{
	return std::make_unique<mpi::stencil_op<xxvii_pt>>(topo);
}


template<>
void copy_services<multilevel_wrapper<solver<xxvii_pt>>>(multilevel_wrapper<solver<xxvii_pt>> & slv, service_manager<stypes> & oldsman)
{
	service_manager<stypes> & newsman = slv.get_services();

	newsman.set_user_reg(oldsman.get_user_reg());
	newsman.set<mpi::message_passing>(oldsman.get_key<mpi::message_passing>());
	newsman.set<mpi::halo_exchange>(oldsman.get_key<mpi::halo_exchange>());
}

}}}
