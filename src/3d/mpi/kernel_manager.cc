#include <cedar/3d/mpi/relax.h>
#include <cedar/3d/mpi/relax_planes.h>
#include <cedar/3d/mpi/coarsen.h>
#include <cedar/3d/mpi/interp.h>
#include <cedar/3d/mpi/residual.h>
#include <cedar/3d/mpi/restrict.h>
#include <cedar/3d/mpi/matvec.h>
#include <cedar/3d/mpi/setup_nog.h>
#include <cedar/3d/mpi/msg_exchanger.h>

#include <cedar/3d/mpi/kernel_manager.h>

namespace cedar { namespace cdr3 { namespace mpi {

kman_ptr build_kernel_manager(config & conf)
{
	return build_kernel_manager(build_kernel_params(conf));
}


kman_ptr build_kernel_manager(std::shared_ptr<kernel_params> params)
{
	auto kman = std::make_unique<kernel_manager<klist<stypes, exec_mode::mpi>>>(params);
	kman->add<point_relax, rbgs>("system");
	kman->add<plane_relax<relax_dir::xy>, planes<relax_dir::xy>>("system");
	kman->add<plane_relax<relax_dir::xz>, planes<relax_dir::xz>>("system");
	kman->add<plane_relax<relax_dir::yz>, planes<relax_dir::yz>>("system");
	kman->add<coarsen_op, galerkin>("system");
	kman->add<interp_add, interp_f90>("system");
	kman->add<restriction, restrict_f90>("system");
	kman->add<residual, residual_f90>("system");
	kman->add<setup_interp, setup_interp_f90>("system");
	kman->add<halo_exchange, msg_exchanger>("msg");
	kman->add<setup_nog, setup_nog_f90>("system");
	kman->add<matvec, matvec_f90>("system");

	kman->set<point_relax>("system");
	kman->set<plane_relax<relax_dir::xy>>("system");
	kman->set<plane_relax<relax_dir::xz>>("system");
	kman->set<plane_relax<relax_dir::yz>>("system");
	kman->set<coarsen_op>("system");
	kman->set<interp_add>("system");
	kman->set<restriction>("system");
	kman->set<residual>("system");
	kman->set<setup_interp>("system");
	kman->set<halo_exchange>("msg");
	kman->set<setup_nog>("system");
	kman->set<matvec>("system");

	// run this manually for now
	auto hex = kman->get_ptr<halo_exchange>();
	kman->get<point_relax>().add_halo(hex.get());
	kman->get<plane_relax<relax_dir::xy>>().add_halo(hex.get());
	kman->get<plane_relax<relax_dir::xz>>().add_halo(hex.get());
	kman->get<plane_relax<relax_dir::yz>>().add_halo(hex.get());
	kman->get<coarsen_op>().add_halo(hex.get());
	kman->get<interp_add>().add_halo(hex.get());
	kman->get<restriction>().add_halo(hex.get());
	kman->get<residual>().add_halo(hex.get());
	kman->get<setup_interp>().add_halo(hex.get());
	kman->get<matvec>().add_halo(hex.get());
	kman->get<setup_nog>().add_halo(hex.get());

	return kman;
}

}}}
