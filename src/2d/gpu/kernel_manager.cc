#include <cedar/2d/gpu/relax.h>
#include <cedar/2d/gpu/coarsen.h>
#include <cedar/2d/gpu/interp.h>
#include <cedar/2d/gpu/residual.h>
#include <cedar/2d/gpu/restrict.h>
#include <cedar/2d/gpu/matvec.h>
#include <cedar/2d/gpu/setup_nog.h>
#include <cedar/2d/gpu/solve_cg.h>
#include <cedar/2d/gpu/msg_exchanger.h>
#ifdef WITH_TAUSCH
#include <cedar/2d/gpu/tausch_exchanger.h>
#endif
#include <cedar/mpi/mpi_wrapper.h>
#include <cedar/malloc_pool.h>

#include <cedar/2d/gpu/kernel_manager.h>

namespace cedar { namespace cdr2 { namespace gpu { namespace mpi {

kman_ptr build_kernel_manager(config & conf)
{
	return build_kernel_manager(build_kernel_params(conf));
}


kman_ptr build_kernel_manager(std::shared_ptr<kernel_params> params)
{
	log::status << *params << std::endl;

	auto kman = std::make_shared<kernel_manager<klist<stypes, exec_mode::mpi>,
	                                            stypes>>(params);
	kman->add<point_relax, rbgs>("system");
	kman->add<coarsen_op, galerkin>("system");
	kman->add<interp_add, interp_f90>("system");
	kman->add<restriction, restrict_f90>("system");
	kman->add<residual, residual_f90>("system");
	kman->add<setup_interp, setup_interp_f90>("system");
	kman->add<setup_nog, setup_nog_f90>("system");
	kman->add<matvec, matvec_f90>("system");
        kman->add<solve_cg, solve_cg_f90>("system");

	kman->set<point_relax>("system");
	kman->set<coarsen_op>("system");
	kman->set<interp_add>("system");
	kman->set<restriction>("system");
	kman->set<residual>("system");
	kman->set<setup_interp>("system");
	kman->set<setup_nog>("system");
	kman->set<matvec>("system");
        kman->set<solve_cg>("system");

	// register services
	auto & services = kman->services();
	services.add<mempool, malloc_pool>("malloc");
	services.add<message_passing, mpi_wrapper>("mpi");
	services.add<halo_exchange, msg_exchanger>("msg");
	#ifdef WITH_TAUSCH
	services.add<halo_exchange, tausch_exchanger>("tausch");
	#endif
	std::string halo_name("msg");
	if (params->halo == kernel_params::halo_lib::tausch)
		halo_name = "tausch";
	services.set<halo_exchange>(halo_name);
	services.set<message_passing>("mpi");
	services.set<mempool>("malloc");

	return kman;
}

}}}}
