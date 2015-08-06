#include "name.h"
#include "residual.h"
#include "setup_interp.h"
#include "setup_relax.h"
#include "galerkin_prod.h"
#include "setup_cg_lu.h"
#include "solve_cg.h"
#include "relax.h"
#include "interp.h"
#include "restrict.h"
#include "setup_nog.h"
#include "halo.h"

#include "manager.h"

using namespace boxmg::bmg2d::kernel;

std::unique_ptr<boxmg::KernelManager> Manager::instance = nullptr;
std::unique_ptr<std::map<std::string,boxmg::KernelManager>> Manager::avail = nullptr;


void Manager::init_reg()
{
	avail = std::make_unique<std::map<std::string, boxmg::KernelManager>>();

	avail->emplace(std::make_pair(name::residual, boxmg::KernelManager()));
	(*avail)[name::residual].add("fortran",
	                             Kernel<const core::StencilOp &,
	                             const core::GridFunc &,
	                             const core::GridFunc &,
	                             core::GridFunc&>(impls::residual_fortran));

	(*avail)[name::residual].add("fortran-msg",
	                             Kernel<const core::StencilOp &,
	                             const core::GridFunc &,
	                             const core::GridFunc &,
	                             core::GridFunc&>(impls::mpi_residual_fortran));

	(*avail)[name::residual].add("c++",
	                             Kernel<const core::StencilOp &,
	                             const core::GridFunc &,
	                             const core::GridFunc &,
	                             core::GridFunc&>(impls::residual));

	(*avail)[name::setup_interp].add("fortran",
	                                 Kernel<int, int , int,
	                                 const core::StencilOp &,
	                                 const core::StencilOp &,
	                                 inter::ProlongOp &>(impls::setup_interp));
	(*avail)[name::galerkin_prod].add("fortran-ex",
	                                  Kernel<int,int,int,
	                                  const inter::ProlongOp&,
	                                  const core::StencilOp&,
	                                  core::StencilOp&>(impls::galerkin_prod));
	(*avail)[name::setup_relax].add("fortran-rbgs-point",
	                                Kernel<const core::StencilOp&,
	                                core::RelaxStencil&>(impls::setup_rbgs_point));
	(*avail)[name::setup_cg_lu].add("fortran",
	                                Kernel<const core::StencilOp&,
	                                core::GridFunc&>(impls::setup_cg_lu));
	(*avail)[name::relax].add("fortran-rbgs",
	                          Kernel<const core::StencilOp&,
	                          core::GridFunc&,
	                          const core::GridFunc&,
	                          const core::RelaxStencil&,
	                          cycle::Dir>(impls::relax_rbgs_point));
	(*avail)[name::relax].add("fortran-msg-rbgs",
	                          Kernel<const core::StencilOp&,
	                          core::GridFunc&,
	                          const core::GridFunc&,
	                          const core::RelaxStencil&,
	                          cycle::Dir>(impls::mpi_relax_rbgs_point));
	(*avail)[name::restriction].add("fortran",
	                                Kernel<const inter::RestrictOp&,
	                                const core::GridFunc&,
	                                core::GridFunc&>(impls::fortran_restrict));
	(*avail)[name::restriction].add("fortran-msg",
	                                Kernel<const inter::RestrictOp&,
	                                const core::GridFunc&,
	                                core::GridFunc&>(impls::mpi_fortran_restrict));
	(*avail)[name::interp_add].add("fortran",
	                               Kernel<const inter::ProlongOp&,
	                               const core::GridFunc&,
	                               const core::GridFunc&,
	                               core::GridFunc&>(impls::fortran_interp));
	(*avail)[name::interp_add].add("fortran-msg",
	                               Kernel<const inter::ProlongOp&,
	                               const core::GridFunc&,
	                               const core::GridFunc&,
	                               core::GridFunc&>(impls::mpi_fortran_interp));
	(*avail)[name::solve_cg].add("fortran",
	                             Kernel<core::GridFunc&,
	                             const core::GridFunc&,
	                             const core::GridFunc&>(impls::fortran_solve_cg));
	(*avail)[name::setup_nog].add("fortran",
	                              Kernel<core::mpi::GridTopo&,
	                              len_t, int*>(impls::fortran_setup_nog));
	(*avail)[name::halo_setup].add("fortran-msg",
	                               Kernel<core::mpi::GridTopo&,void**>(impls::setup_msg));
	(*avail)[name::halo_exchange].add("fortran-msg",
	                                  Kernel<core::mpi::GridFunc&>(impls::msg_exchange));
	(*avail)[name::halo_stencil_exchange].add("fortran-msg",
	                                          Kernel<core::mpi::StencilOp&>(impls::msg_stencil_exchange));
	(*avail)[name::setup_interp].add("fortran-msg",
	                                 Kernel<int,int,int,
	                                 const core::StencilOp&,
	                                 const core::StencilOp&,
	                                 inter::ProlongOp&>(impls::mpi_setup_interp));
	(*avail)[name::galerkin_prod].add("fortran-msg",
	                                  Kernel<int,int,int,
	                                  const inter::ProlongOp&,
	                                  const core::StencilOp&,
	                                  core::StencilOp&>(impls::mpi_galerkin_prod));
	(*avail)[name::setup_relax].add("fortran-msg-rbgs-point",
	                                Kernel<const core::StencilOp&,
	                                core::RelaxStencil&>(impls::mpi_setup_rbgs_point));
}


void Manager::init()
{
	if (!avail) init_reg();
	instance = std::make_unique<boxmg::KernelManager>();

	using namespace boxmg;
	using namespace boxmg::bmg2d;

	// TODO: This could be automatic.

	std::string res_kern = config::get<std::string>("kernels.residual", "c++");
	std::string setup_interp_kern = config::get<std::string>("kernels.setup-interp", "fortran");
	std::string galerkin_prod_kern = config::get<std::string>("kernels.galerkin-prod", "fortran-ex");
	std::string setup_relax_kern = config::get<std::string>("kernels.setup-relax", "fortran-rbgs-point");
	std::string setup_cg_lu_kern = config::get<std::string>("kernels.setup-cg-lu", "fortran");
	std::string relax_kern = config::get<std::string>("kernels.relax", "fortran-rbgs");
	std::string restrict_kern = config::get<std::string>("kernels.restrict", "fortran");
	std::string interp_kern = config::get<std::string>("kernels.interp-add", "fortran");
	std::string solve_cg_kern = config::get<std::string>("kernels.solve-cg", "fortran");
	std::string halo_setup_kern = config::get<std::string>("kernels.halo-setup", "fortran-msg");
	std::string halo_exchange_kern = config::get<std::string>("kernels.halo-exchange", "fortran-msg");

	log::info << "Using '" << res_kern << "' to compute residual." << std::endl;
	log::info << "Using '" + setup_interp_kern + "' to setup interpolation operator." << std::endl;
	log::info << "Using '" + galerkin_prod_kern + "' for operator coarsening." << std::endl;
	log::info << "Using '" + setup_relax_kern + "' to setup relaxation." << std::endl;
	log::info << "Using '" + setup_cg_lu_kern + "' to setup coarse grid solve." << std::endl;
	log::info << "Using '" + relax_kern + "' for relaxation." << std::endl;
	log::info << "Using '" + restrict_kern + "' for " << name::restriction << "." << std::endl;
	log::info << "Using '" + interp_kern + "' for " << name::interp_add << "." << std::endl;
	log::info << "Using '" + solve_cg_kern + "' for " << name::solve_cg << "." << std::endl;

	(*instance)[name::residual] = (*avail)[name::residual].at(res_kern);
	(*instance)[name::setup_interp] = (*avail)[name::setup_interp].at(setup_interp_kern);
	(*instance)[name::galerkin_prod] = (*avail)[name::galerkin_prod].at(galerkin_prod_kern);
	(*instance)[name::setup_relax] = (*avail)[name::setup_relax].at(setup_relax_kern);
	(*instance)[name::setup_cg_lu] = (*avail)[name::setup_cg_lu].at(setup_cg_lu_kern);
	(*instance)[name::relax] = (*avail)[name::relax].at(relax_kern);
	(*instance)[name::restriction] = (*avail)[name::restriction].at(restrict_kern);
	(*instance)[name::interp_add] = (*avail)[name::interp_add].at(restrict_kern);
	(*instance)[name::solve_cg] = (*avail)[name::solve_cg].at(solve_cg_kern);
	(*instance)[name::setup_nog] = (*avail)[name::setup_nog].at("fortran");
	(*instance)[name::halo_setup] = (*avail)[name::halo_setup].at("fortran-msg");
	(*instance)[name::halo_exchange] = (*avail)[name::halo_exchange].at("fortran-msg");
	(*instance)[name::halo_stencil_exchange] = (*avail)[name::halo_stencil_exchange].at("fortran-msg");
	(*instance)[name::setup_interp] = (*avail)[name::setup_interp].at("fortran-msg");
	(*instance)[name::galerkin_prod] = (*avail)[name::galerkin_prod].at("fortran-msg");
	(*instance)[name::setup_relax] = (*avail)[name::setup_relax].at("fortran-msg-rbgs-point");
	(*instance)[name::relax] = (*avail)[name::relax].at("fortran-msg-rbgs");
	(*instance)[name::restriction] = (*avail)[name::restriction].at("fortran-msg");
	(*instance)[name::interp_add] = (*avail)[name::interp_add].at("fortran-msg");
	(*instance)[name::residual] = (*avail)[name::residual].at("fortran-msg");
}


void Manager::check_init()
{
	if (!instance)
		init();
}




void Manager::change(const std::string & kname, const std::string & rname)
{
	check_init();
	(*instance)[kname] = (*avail)[kname][rname];
	log::info << "Now using '" << rname << "' to compute " << kname << '.' << std::endl;
}
