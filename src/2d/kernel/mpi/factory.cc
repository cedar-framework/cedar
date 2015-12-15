#include <memory>
#include "boxmg/2d/kernel/name.h"
#include "boxmg/2d/residual.h"
#include "boxmg/2d/inter/setup_interp.h"
#include "boxmg/2d/relax/setup_relax.h"
#include "boxmg/2d/inter/galerkin_prod.h"
#include "boxmg/2d/cg/setup_cg_lu.h"
#include "boxmg/2d/cg/solve_cg.h"
#include "boxmg/2d/relax/relax.h"
#include "boxmg/2d/inter/interp.h"
#include "boxmg/2d/inter/restrict.h"
#include "boxmg/2d/kernel/setup_nog.h"
#include "boxmg/2d/cg/setup_cg_boxmg.h"
#include "boxmg/2d/mpi/halo.h"
#include "boxmg/2d/matvec.h"

#include "boxmg/2d/kernel/mpi/factory.h"

namespace boxmg { namespace bmg2d {
			class BoxMG;
}}

namespace boxmg { namespace bmg2d { namespace kernel { namespace mpi {


namespace factory
{
	std::shared_ptr<Registry> from_config(config::Reader &conf)
	{
		auto kreg = std::make_shared<Registry>();

		kreg->add(name::residual, "fortran-msg",
		         Kernel<const bmg2d::stencil_op &,
		         const bmg2d::grid_func &,
		         const bmg2d::grid_func &,
		         bmg2d::grid_func&>(impls::mpi_residual_fortran));

		kreg->add(name::setup_cg_lu, "fortran-msg",
		         Kernel<const bmg2d::stencil_op&,
		         bmg2d::grid_func&>(impls::mpi_setup_cg_lu));

		kreg->add(name::relax, "fortran-msg-rbgs",
		         Kernel<const bmg2d::stencil_op&,
		         bmg2d::grid_func&,
		         const bmg2d::grid_func&,
		         const bmg2d::relax_stencil&,
		         cycle::Dir>(impls::mpi_relax_rbgs_point));

		kreg->add(name::relax_lines_x, "fortran-msg",
		         Kernel<const bmg2d::stencil_op&,
		         bmg2d::grid_func&,
		         const bmg2d::grid_func&,
		         const bmg2d::relax_stencil&,
		          bmg2d::grid_func&,
		         cycle::Dir>(impls::mpi_relax_lines_x));

		kreg->add(name::relax_lines_y, "fortran-msg",
		         Kernel<const bmg2d::stencil_op&,
		         bmg2d::grid_func&,
		         const bmg2d::grid_func&,
		         const bmg2d::relax_stencil&,
		          bmg2d::grid_func&,
		         cycle::Dir>(impls::mpi_relax_lines_y));

		kreg->add(name::restriction, "fortran-msg",
		         Kernel<const inter::restrict_op&,
		         const bmg2d::grid_func&,
		         bmg2d::grid_func&>(impls::mpi_fortran_restrict));

		kreg->add(name::interp_add, "fortran-msg",
		         Kernel<const inter::prolong_op&,
		         const bmg2d::grid_func&,
		         const bmg2d::grid_func&,
		         bmg2d::grid_func&>(impls::mpi_fortran_interp));

		kreg->add(name::solve_cg, "fortran-msg",
		          Kernel<bmg2d::grid_func&,
		          const bmg2d::grid_func&,
		          const bmg2d::grid_func&,
		          real_t*>(impls::mpi_solve_cg_lu));

		kreg->add(name::setup_nog, "fortran",
		         Kernel<bmg2d::mpi::grid_topo&,
		         len_t, int*>(impls::fortran_setup_nog));

		kreg->add(name::halo_setup, "fortran-msg",
		         Kernel<bmg2d::mpi::grid_topo&,void**>(impls::setup_msg));

		kreg->add(name::halo_exchange, "fortran-msg",
		         Kernel<bmg2d::mpi::grid_func&>(impls::msg_exchange));

		kreg->add(name::halo_stencil_exchange, "fortran-msg",
		         Kernel<bmg2d::mpi::stencil_op&>(impls::msg_stencil_exchange));

		kreg->add(name::setup_interp, "fortran-msg",
		         Kernel<int,int,int,
		         const bmg2d::stencil_op&,
		         const bmg2d::stencil_op&,
		         inter::prolong_op&>(impls::mpi_setup_interp));

		kreg->add(name::galerkin_prod, "fortran-msg",
		         Kernel<int,int,int,
		         const inter::prolong_op&,
		         const bmg2d::stencil_op&,
		         bmg2d::stencil_op&>(impls::mpi_galerkin_prod));

		kreg->add(name::setup_relax, "fortran-msg-rbgs-point",
		         Kernel<const bmg2d::stencil_op&,
		         bmg2d::relax_stencil&>(impls::mpi_setup_rbgs_point));

		kreg->add(name::setup_relax_x, "fortran-msg",
		         Kernel<const bmg2d::stencil_op&,
		         bmg2d::relax_stencil&>(impls::mpi_setup_rbgs_x));

		kreg->add(name::setup_relax_y, "fortran-msg",
		         Kernel<const bmg2d::stencil_op&,
		         bmg2d::relax_stencil&>(impls::mpi_setup_rbgs_y));

		kreg->add(name::setup_cg_boxmg, "fortran-msg",
		          Kernel<const bmg2d::stencil_op &,
		          std::shared_ptr<solver>*>(impls::setup_cg_boxmg));

		kreg->add(name::solve_cg_boxmg, "fortran-msg",
		          Kernel<const solver &,
		          bmg2d::grid_func &,
		          const bmg2d::grid_func &>(impls::solve_cg_boxmg));

		kreg->add(name::matvec, "fortran-msg",
		          Kernel<const bmg2d::stencil_op&,
		          const bmg2d::grid_func &, bmg2d::grid_func&>(impls::matvec));


		std::vector<std::tuple<std::string, std::string, std::string>> defaults = {
			std::make_tuple(name::residual, "kernels.residual", "fortran-msg"),
			std::make_tuple(name::setup_interp, "kernels.setup-interp", "fortran-msg"),
			std::make_tuple(name::galerkin_prod, "kernels.galerkin-prod", "fortran-msg"),
			std::make_tuple(name::setup_relax, "kernels.setup-relax", "fortran-msg-rbgs-point"),
			std::make_tuple(name::setup_cg_lu, "kernels.setup-cg-lu", "fortran-msg"),
			std::make_tuple(name::restriction, "kernels.restrict", "fortran-msg"),
			std::make_tuple(name::interp_add, "kernels.interp-add", "fortran-msg"),
			std::make_tuple(name::solve_cg, "kernels.solve-cg", "fortran-msg"),
			std::make_tuple(name::relax, "kernels.relax", "fortran-msg-rbgs"),
			std::make_tuple(name::setup_relax_x, "kernels.setup-relax-x", "fortran-msg"),
			std::make_tuple(name::setup_relax_y, "kernels.setup-relax-y", "fortran-msg"),
			std::make_tuple(name::relax_lines_x, "kernels.relax-x", "fortran-msg"),
			std::make_tuple(name::relax_lines_y, "kernels.relax-y", "fortran-msg"),
			std::make_tuple(name::setup_nog, "kernels.setup-nog", "fortran"),
			std::make_tuple(name::halo_setup, "kernels.halo-setup", "fortran-msg"),
			std::make_tuple(name::halo_exchange, "kernels.halo-exchange", "fortran-msg"),
			std::make_tuple(name::halo_stencil_exchange, "kernels.halo-stencil-exchange", "fortran-msg"),
			std::make_tuple(name::matvec, "kernels.matvec", "fortran-msg"),
			std::make_tuple(name::setup_cg_boxmg, "kernels.setup-cg-boxmg", "fortran-msg"),
			std::make_tuple(name::solve_cg_boxmg, "kernels.solve-cg-boxmg", "fortran-msg")
		};

		for (auto&& v : defaults) {
			std::string kname = conf.get<std::string>(std::get<1>(v), std::get<2>(v));
			log::info << "Using '" + kname + " ' for " <<  std::get<0>(v) << "." << std::endl;
			kreg->set(std::get<0>(v), kname);
		}

		return kreg;
	}
}

}}}}
