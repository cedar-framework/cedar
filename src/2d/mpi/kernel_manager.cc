#include <cedar/2d/mpi/relax.h>
#include <cedar/2d/mpi/ml_relax.h>
#include <cedar/2d/mpi/coarsen.h>
#include <cedar/2d/mpi/interp.h>
#include <cedar/2d/mpi/residual.h>
#include <cedar/2d/mpi/restrict.h>
#include <cedar/2d/mpi/matvec.h>
#include <cedar/2d/mpi/setup_nog.h>
#include <cedar/2d/mpi/msg_exchanger.h>
#ifdef WITH_TAUSCH
#include <cedar/2d/mpi/tausch_exchanger.h>
#endif
#include <cedar/mpi/mpi_wrapper.h>
#include <cedar/malloc_pool.h>

#include <cedar/2d/mpi/kernel_manager.h>

namespace cedar { namespace cdr2 { namespace mpi {

kman_ptr build_kernel_manager(config & conf)
{
	return build_kernel_manager(build_kernel_params(conf));
}


kman_ptr build_kernel_manager(std::shared_ptr<kernel_params> params)
{
	log::status << *params << std::endl;

	auto kman = std::make_shared<kernel_manager<klist<stypes, exec_mode::mpi>,
	                                            stypes>>(params);
        /* CPU kernels */
	kman->add<point_relax, rbgs<cedar::cpu>>("system");
	kman->add<line_relax<relax_dir::x>, lines<relax_dir::x>>("two-level");
	kman->add<line_relax<relax_dir::y>, lines<relax_dir::y>>("two-level");
	kman->add<line_relax<relax_dir::x>, ml_line_relax<relax_dir::x>>("n-level");
	kman->add<line_relax<relax_dir::y>, ml_line_relax<relax_dir::y>>("n-level");
	kman->add<line_relax<relax_dir::x>, ml_line_relax<relax_dir::x>>("n-level-elim", false);
	kman->add<line_relax<relax_dir::y>, ml_line_relax<relax_dir::y>>("n-level-elim", false);
	kman->add<coarsen_op, galerkin<cedar::cpu>>("system");
	kman->add<interp_add, interp_f90<cedar::cpu>>("system");
	kman->add<restriction, restrict_f90<cedar::cpu>>("system");
	kman->add<residual, residual_f90<cedar::cpu>>("system");
	kman->add<setup_interp, setup_interp_f90<cedar::cpu>>("system");
	kman->add<setup_nog, setup_nog_f90>("system");
	kman->add<matvec, matvec_f90<cedar::cpu>>("system");

        /* GPU kernels */
        kman->add<point_relax, rbgs<cedar::gpu>>("gpu");
	kman->add<coarsen_op, galerkin<cedar::gpu>>("gpu");
	kman->add<interp_add, interp_f90<cedar::gpu>>("gpu");
	kman->add<restriction, restrict_f90<cedar::gpu>>("gpu");
	kman->add<residual, residual_f90<cedar::gpu>>("gpu");
	kman->add<setup_interp, setup_interp_f90<cedar::gpu>>("gpu");
	kman->add<matvec, matvec_f90<cedar::gpu>>("gpu");

	std::string line_relax_name("two-level");
	if (params->ml_relax.enabled) {
		if (params->ml_relax.factorize)
			line_relax_name = "n-level";
		else
			line_relax_name = "n-level-elim";
	}
	log::debug << "Using <" << line_relax_name << "> for line relaxation" << std::endl;
	kman->set<line_relax<relax_dir::x>>(line_relax_name);
	kman->set<line_relax<relax_dir::y>>(line_relax_name);

        /* Select between CPU and GPU kernels */
        const std::string kernels = (params->use_gpu ? "gpu" : "system");
        kman->set<point_relax>(kernels);
	kman->set<coarsen_op>(kernels);
	kman->set<interp_add>(kernels);
	kman->set<restriction>(kernels);
	kman->set<residual>(kernels);
	kman->set<setup_interp>(kernels);
	kman->set<setup_nog>("system");
	kman->set<matvec>(kernels);

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

}}}
