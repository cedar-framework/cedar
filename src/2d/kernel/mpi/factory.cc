#include <memory>
#include "cedar/kernel.h"
#include "cedar/kernel_name.h"
#include "cedar/2d/residual.h"
#include "cedar/2d/inter/setup_interp.h"
#include "cedar/2d/relax/setup_relax.h"
#include "cedar/2d/inter/galerkin_prod.h"
#include "cedar/2d/cg/setup_cg_lu.h"
#include "cedar/2d/cg/solve_cg.h"
#include "cedar/2d/relax/relax.h"
#include "cedar/2d/inter/interp.h"
#include "cedar/2d/inter/restrict.h"
#include "cedar/2d/kernel/setup_nog.h"
#include "cedar/2d/cg/setup_cg_boxmg.h"
#include "cedar/2d/cg/setup_cg_redist.h"
#include "cedar/2d/mpi/halo.h"
#include "cedar/2d/matvec.h"

#include "cedar/2d/kernel/mpi/factory.h"

namespace cedar { namespace cdr2 {
			class BoxMG;
}}

namespace cedar { namespace cdr2 { namespace kernel { namespace mpi {

namespace mpi = cedar::cdr2::mpi;

namespace factory
{
	void init(std::shared_ptr<registry> kreg, config::reader & conf)
	{
		auto params = build_kernel_params(conf);

		kreg->add(kernel_name::residual, "fortran-msg",
		         cedar::kernel<const mpi::stencil_op &,
		          const mpi::grid_func &,
		          const mpi::grid_func &,
		          mpi::grid_func&>(impls::mpi_residual_fortran, params));

		kreg->add(kernel_name::setup_cg_lu, "fortran-msg",
		         cedar::kernel<const mpi::stencil_op&,
		         mpi::grid_func&>(impls::mpi_setup_cg_lu, params));

		kreg->add(kernel_name::relax, "fortran-msg-rbgs",
		          cedar::kernel<const mpi::stencil_op&,
		          mpi::grid_func&,
		          const mpi::grid_func&,
		          const cdr2::relax_stencil&,
		          cycle::Dir>(impls::mpi_relax_rbgs_point, params));

		kreg->add(kernel_name::relax_lines_x, "fortran-msg",
		          cedar::kernel<const mpi::stencil_op&,
		          mpi::grid_func&,
		          const mpi::grid_func&,
		          const cdr2::relax_stencil&,
		          mpi::grid_func&,
		          cycle::Dir>(impls::mpi_relax_lines_x, params));

		kreg->add(kernel_name::relax_lines_y, "fortran-msg",
		          cedar::kernel<const mpi::stencil_op&,
		          mpi::grid_func&,
		          const mpi::grid_func&,
		          const cdr2::relax_stencil&,
		          mpi::grid_func&,
		         cycle::Dir>(impls::mpi_relax_lines_y, params));

		kreg->add(kernel_name::restriction, "fortran-msg",
		          cedar::kernel<const inter::mpi::restrict_op&,
		          const mpi::grid_func&,
		          mpi::grid_func&>(impls::mpi_fortran_restrict, params));

		kreg->add(kernel_name::interp_add, "fortran-msg",
		          cedar::kernel<const inter::mpi::prolong_op&,
		         const mpi::grid_func&,
		         const mpi::grid_func&,
		         mpi::grid_func&>(impls::mpi_fortran_interp, params));

		kreg->add(kernel_name::solve_cg, "fortran-msg",
		          cedar::kernel<mpi::grid_func&,
		          const mpi::grid_func&,
		          const mpi::grid_func&,
		          real_t*>(impls::mpi_solve_cg_lu, params));

		kreg->add(kernel_name::setup_nog, "fortran",
		         cedar::kernel<grid_topo&,
		         len_t, int*>(impls::fortran_setup_nog, params));

		kreg->add(kernel_name::halo_setup, "fortran-msg",
		         cedar::kernel<grid_topo&,void**>(impls::setup_msg, params));

		kreg->add(kernel_name::halo_exchange, "fortran-msg",
		          cedar::kernel<cdr2::mpi::grid_func&>(impls::msg_exchange, params));

		kreg->add(kernel_name::halo_stencil_exchange, "fortran-msg",
		          cedar::kernel<cdr2::mpi::stencil_op&>(impls::msg_stencil_exchange, params));

		kreg->add(kernel_name::setup_interp, "fortran-msg",
		          cedar::kernel<int,int,int,
		          const mpi::stencil_op&,
		          const mpi::stencil_op&,
		          inter::mpi::prolong_op&>(impls::mpi_setup_interp, params));

		kreg->add(kernel_name::galerkin_prod, "fortran-msg",
		          cedar::kernel<int,int,int,
		          const inter::mpi::prolong_op&,
		          const mpi::stencil_op&,
		          mpi::stencil_op&>(impls::mpi_galerkin_prod, params));

		kreg->add(kernel_name::setup_relax, "fortran-msg-rbgs-point",
		         cedar::kernel<const mpi::stencil_op&,
		         cdr2::relax_stencil&>(impls::mpi_setup_rbgs_point, params));

		kreg->add(kernel_name::setup_relax_x, "fortran-msg",
		          cedar::kernel<const mpi::stencil_op&,
		          cdr2::relax_stencil&>(impls::mpi_setup_rbgs_x, params));

		kreg->add(kernel_name::setup_relax_y, "fortran-msg",
		         cedar::kernel<const mpi::stencil_op&,
		         cdr2::relax_stencil&>(impls::mpi_setup_rbgs_y, params));

		kreg->add(kernel_name::setup_cg_boxmg, "fortran-msg",
		          cedar::kernel<const mpi::stencil_op &,
		          std::shared_ptr<config::reader>,
		          std::shared_ptr<solver>*>(impls::setup_cg_boxmg, params));

		kreg->add(kernel_name::solve_cg_boxmg, "fortran-msg",
		          cedar::kernel<const solver &,
		          mpi::grid_func &,
		          const mpi::grid_func &>(impls::solve_cg_boxmg, params));

		kreg->add(kernel_name::matvec, "fortran-msg",
		          cedar::kernel<const mpi::stencil_op&,
		          const mpi::grid_func &, mpi::grid_func&>(impls::matvec, params));

		kreg->add(kernel_name::setup_cg_redist, "c++",
		          cedar::kernel<const mpi::stencil_op &,
		          std::shared_ptr<config::reader>,
		          std::shared_ptr<mpi::redist_solver>*,
		          std::vector<int>&>(impls::setup_cg_redist, params));

		kreg->add(kernel_name::solve_cg_redist, "c++",
		          cedar::kernel<const mpi::redist_solver &,
		          mpi::grid_func &,
		          const mpi::grid_func &>(impls::solve_cg_redist, params));


		std::vector<std::tuple<std::string, std::string, std::string>> defaults = {
			std::make_tuple(kernel_name::residual, "kernels.residual", "fortran-msg"),
			std::make_tuple(kernel_name::setup_interp, "kernels.setup-interp", "fortran-msg"),
			std::make_tuple(kernel_name::galerkin_prod, "kernels.galerkin-prod", "fortran-msg"),
			std::make_tuple(kernel_name::setup_relax, "kernels.setup-relax", "fortran-msg-rbgs-point"),
			std::make_tuple(kernel_name::setup_cg_lu, "kernels.setup-cg-lu", "fortran-msg"),
			std::make_tuple(kernel_name::restriction, "kernels.restrict", "fortran-msg"),
			std::make_tuple(kernel_name::interp_add, "kernels.interp-add", "fortran-msg"),
			std::make_tuple(kernel_name::solve_cg, "kernels.solve-cg", "fortran-msg"),
			std::make_tuple(kernel_name::relax, "kernels.relax", "fortran-msg-rbgs"),
			std::make_tuple(kernel_name::setup_relax_x, "kernels.setup-relax-x", "fortran-msg"),
			std::make_tuple(kernel_name::setup_relax_y, "kernels.setup-relax-y", "fortran-msg"),
			std::make_tuple(kernel_name::relax_lines_x, "kernels.relax-x", "fortran-msg"),
			std::make_tuple(kernel_name::relax_lines_y, "kernels.relax-y", "fortran-msg"),
			std::make_tuple(kernel_name::setup_nog, "kernels.setup-nog", "fortran"),
			std::make_tuple(kernel_name::halo_setup, "kernels.halo-setup", "fortran-msg"),
			std::make_tuple(kernel_name::halo_exchange, "kernels.halo-exchange", "fortran-msg"),
			std::make_tuple(kernel_name::halo_stencil_exchange, "kernels.halo-stencil-exchange", "fortran-msg"),
			std::make_tuple(kernel_name::matvec, "kernels.matvec", "fortran-msg"),
			std::make_tuple(kernel_name::setup_cg_boxmg, "kernels.setup-cg-boxmg", "fortran-msg"),
			std::make_tuple(kernel_name::solve_cg_boxmg, "kernels.solve-cg-boxmg", "fortran-msg"),
			std::make_tuple(kernel_name::setup_cg_redist, "kernels.setup-cg-redist", "c++"),
			std::make_tuple(kernel_name::solve_cg_redist, "kernels.solve-cg-redist", "c++")
		};

		for (auto&& v : defaults) {
			std::string kname = conf.get<std::string>(std::get<1>(v), std::get<2>(v));
			log::debug << "Using '" + kname + "' for " <<  std::get<0>(v) << "." << std::endl;
			kreg->set(std::get<0>(v), kname);
		}
	}


	std::shared_ptr<registry> from_config(config::reader &conf)
	{
		auto kreg = std::make_shared<registry>();

		init(kreg, conf);

		return kreg;
	}
}

}}}}
