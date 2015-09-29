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

#include "factory.h"

namespace boxmg { namespace bmg2d { namespace kernel {


namespace factory
{
	std::shared_ptr<Registry> from_config(config::Reader &conf)
	{
		auto kreg = std::make_shared<Registry>();


		kreg->add(name::residual, "fortran",
		         Kernel<const core::StencilOp &,
		         const core::GridFunc &,
		         const core::GridFunc &,
		         core::GridFunc&>(impls::residual_fortran));

		kreg->add(name::residual, "fortran-msg",
		         Kernel<const core::StencilOp &,
		         const core::GridFunc &,
		         const core::GridFunc &,
		         core::GridFunc&>(impls::mpi_residual_fortran));

		kreg->add(name::residual,"c++",
		         Kernel<const core::StencilOp &,
		         const core::GridFunc &,
		         const core::GridFunc &,
		         core::GridFunc&>(impls::residual));

		kreg->add(name::setup_interp, "fortran",
		         Kernel<int, int , int,
		         const core::StencilOp &,
		         const core::StencilOp &,
		         inter::ProlongOp &>(impls::setup_interp));

		kreg->add(name::galerkin_prod, "fortran-ex",
		         Kernel<int,int,int,
		         const inter::ProlongOp&,
		         const core::StencilOp&,
		         core::StencilOp&>(impls::galerkin_prod));

		kreg->add(name::setup_relax,"fortran-rbgs-point",
		         Kernel<const core::StencilOp&,
		         core::RelaxStencil&>(impls::setup_rbgs_point));

		kreg->add(name::setup_relax_x,"fortran",
		         Kernel<const core::StencilOp&,
		         core::RelaxStencil&>(impls::setup_rbgs_x));

		kreg->add(name::setup_relax_y,"fortran",
		         Kernel<const core::StencilOp&,
		         core::RelaxStencil&>(impls::setup_rbgs_y));

		kreg->add(name::setup_cg_lu, "fortran",
		         Kernel<const core::StencilOp&,
		         core::GridFunc&>(impls::setup_cg_lu));

		kreg->add(name::relax, "fortran-rbgs",
		         Kernel<const core::StencilOp&,
		         core::GridFunc&,
		         const core::GridFunc&,
		         const core::RelaxStencil&,
		         cycle::Dir>(impls::relax_rbgs_point));

		kreg->add(name::relax_lines_x, "fortran",
		         Kernel<const core::StencilOp&,
		         core::GridFunc&,
		         const core::GridFunc&,
		         const core::RelaxStencil&,
		          core::GridFunc&,
		         cycle::Dir>(impls::relax_lines_x));

		kreg->add(name::relax_lines_y, "fortran",
		         Kernel<const core::StencilOp&,
		         core::GridFunc&,
		         const core::GridFunc&,
		         const core::RelaxStencil&,
		          core::GridFunc&,
		         cycle::Dir>(impls::relax_lines_y));

		kreg->add(name::relax, "fortran-msg-rbgs",
		         Kernel<const core::StencilOp&,
		         core::GridFunc&,
		         const core::GridFunc&,
		         const core::RelaxStencil&,
		         cycle::Dir>(impls::mpi_relax_rbgs_point));

		kreg->add(name::restriction, "fortran",
		         Kernel<const inter::RestrictOp&,
		         const core::GridFunc&,
		         core::GridFunc&>(impls::fortran_restrict));

		kreg->add(name::restriction, "fortran-msg",
		         Kernel<const inter::RestrictOp&,
		         const core::GridFunc&,
		         core::GridFunc&>(impls::mpi_fortran_restrict));

		kreg->add(name::interp_add,"fortran",
		         Kernel<const inter::ProlongOp&,
		         const core::GridFunc&,
		         const core::GridFunc&,
		         core::GridFunc&>(impls::fortran_interp));

		kreg->add(name::interp_add, "fortran-msg",
		         Kernel<const inter::ProlongOp&,
		         const core::GridFunc&,
		         const core::GridFunc&,
		         core::GridFunc&>(impls::mpi_fortran_interp));

		kreg->add(name::solve_cg, "fortran",
		          Kernel<core::GridFunc&,
		          const core::GridFunc&,
		          const core::GridFunc&,
		          real_t*>(impls::fortran_solve_cg));

		kreg->add(name::setup_nog, "fortran",
		         Kernel<core::mpi::GridTopo&,
		         len_t, int*>(impls::fortran_setup_nog));

		kreg->add(name::halo_setup, "fortran-msg",
		         Kernel<core::mpi::GridTopo&,void**>(impls::setup_msg));

		kreg->add(name::halo_exchange, "fortran-msg",
		         Kernel<core::mpi::GridFunc&>(impls::msg_exchange));

		kreg->add(name::halo_stencil_exchange, "fortran-msg",
		         Kernel<core::mpi::StencilOp&>(impls::msg_stencil_exchange));

		kreg->add(name::setup_interp, "fortran-msg",
		         Kernel<int,int,int,
		         const core::StencilOp&,
		         const core::StencilOp&,
		         inter::ProlongOp&>(impls::mpi_setup_interp));

		kreg->add(name::galerkin_prod, "fortran-msg",
		         Kernel<int,int,int,
		         const inter::ProlongOp&,
		         const core::StencilOp&,
		         core::StencilOp&>(impls::mpi_galerkin_prod));

		kreg->add(name::setup_relax, "fortran-msg-rbgs-point",
		         Kernel<const core::StencilOp&,
		         core::RelaxStencil&>(impls::mpi_setup_rbgs_point));

		// TODO: This could be automatic.

		std::string res_kern = conf.get<std::string>("kernels.residual", "c++");
		std::string setup_interp_kern = conf.get<std::string>("kernels.setup-interp", "fortran");
		std::string galerkin_prod_kern = conf.get<std::string>("kernels.galerkin-prod", "fortran-ex");
		std::string setup_relax_kern = conf.get<std::string>("kernels.setup-relax", "fortran-rbgs-point");
		std::string setup_cg_lu_kern = conf.get<std::string>("kernels.setup-cg-lu", "fortran");
		std::string relax_kern = conf.get<std::string>("kernels.relax", "fortran-rbgs");
		std::string restrict_kern = conf.get<std::string>("kernels.restrict", "fortran");
		std::string interp_kern = conf.get<std::string>("kernels.interp-add", "fortran");
		std::string solve_cg_kern = conf.get<std::string>("kernels.solve-cg", "fortran");
		std::string halo_setup_kern = conf.get<std::string>("kernels.halo-setup", "fortran-msg");
		std::string halo_exchange_kern = conf.get<std::string>("kernels.halo-exchange", "fortran-msg");

		log::info << "Using '" << res_kern << "' to compute residual." << std::endl;
		log::info << "Using '" + setup_interp_kern + "' to setup interpolation operator." << std::endl;
		log::info << "Using '" + galerkin_prod_kern + "' for operator coarsening." << std::endl;
		log::info << "Using '" + setup_relax_kern + "' to setup relaxation." << std::endl;
		log::info << "Using '" + setup_cg_lu_kern + "' to setup coarse grid solve." << std::endl;
		log::info << "Using '" + relax_kern + "' for relaxation." << std::endl;
		log::info << "Using '" + restrict_kern + "' for " << name::restriction << "." << std::endl;
		log::info << "Using '" + interp_kern + "' for " << name::interp_add << "." << std::endl;
		log::info << "Using '" + solve_cg_kern + "' for " << name::solve_cg << "." << std::endl;

		kreg->set(name::residual, res_kern);
		kreg->set(name::galerkin_prod, galerkin_prod_kern);
		kreg->set(name::setup_interp, setup_interp_kern);
		kreg->set(name::setup_relax, setup_relax_kern);
		kreg->set(name::setup_relax_x, "fortran");
		kreg->set(name::setup_relax_y, "fortran");
		kreg->set(name::setup_cg_lu, setup_cg_lu_kern);
		kreg->set(name::relax, relax_kern);
		kreg->set(name::restriction, restrict_kern);
		kreg->set(name::interp_add, interp_kern);
		kreg->set(name::solve_cg, solve_cg_kern);
		kreg->set(name::relax_lines_x, "fortran");
		kreg->set(name::relax_lines_y, "fortran");

		return kreg;
	}
}

}}}
