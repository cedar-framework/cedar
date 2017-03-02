#include <memory>
#include "boxmg/kernel.h"
#include "boxmg/kernel_name.h"
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
#include "boxmg/2d/cg/setup_cg_redist.h"
#include "boxmg/2d/mpi/halo.h"
#include "boxmg/2d/matvec.h"

#include "boxmg/2d/kernel/mpi/factory.h"

namespace boxmg { namespace bmg2d {
			class BoxMG;
}}

namespace boxmg { namespace bmg2d { namespace kernel { namespace mpi {

namespace mpi = boxmg::bmg2d::mpi;

namespace factory
{
	std::shared_ptr<registry> from_config(config::reader &conf)
	{
		auto kreg = std::make_shared<registry>();

		kreg->add(kernel_name::residual, "fortran-msg",
		         boxmg::kernel<const mpi::stencil_op &,
		          const mpi::grid_func &,
		          const mpi::grid_func &,
		          mpi::grid_func&>(impls::mpi_residual_fortran));

		kreg->add(kernel_name::setup_cg_lu, "fortran-msg",
		         boxmg::kernel<const mpi::stencil_op&,
		         mpi::grid_func&>(impls::mpi_setup_cg_lu));

		kreg->add(kernel_name::relax, "fortran-msg-rbgs",
		          boxmg::kernel<const mpi::stencil_op&,
		          mpi::grid_func&,
		          const mpi::grid_func&,
		          const bmg2d::relax_stencil&,
		          cycle::Dir>(impls::mpi_relax_rbgs_point));

		kreg->add(kernel_name::relax_lines_x, "fortran-msg",
		          boxmg::kernel<const mpi::stencil_op&,
		          mpi::grid_func&,
		          const mpi::grid_func&,
		          const bmg2d::relax_stencil&,
		          mpi::grid_func&,
		          cycle::Dir>(impls::mpi_relax_lines_x));

		kreg->add(kernel_name::relax_lines_y, "fortran-msg",
		          boxmg::kernel<const mpi::stencil_op&,
		          mpi::grid_func&,
		          const mpi::grid_func&,
		          const bmg2d::relax_stencil&,
		          mpi::grid_func&,
		         cycle::Dir>(impls::mpi_relax_lines_y));

		kreg->add(kernel_name::restriction, "fortran-msg",
		          boxmg::kernel<const inter::mpi::restrict_op&,
		          const mpi::grid_func&,
		          mpi::grid_func&>(impls::mpi_fortran_restrict));

		kreg->add(kernel_name::interp_add, "fortran-msg",
		          boxmg::kernel<const inter::mpi::prolong_op&,
		         const mpi::grid_func&,
		         const mpi::grid_func&,
		         mpi::grid_func&>(impls::mpi_fortran_interp));

		kreg->add(kernel_name::solve_cg, "fortran-msg",
		          boxmg::kernel<mpi::grid_func&,
		          const mpi::grid_func&,
		          const mpi::grid_func&,
		          real_t*>(impls::mpi_solve_cg_lu));

		kreg->add(kernel_name::setup_nog, "fortran",
		         boxmg::kernel<grid_topo&,
		         len_t, int*>(impls::fortran_setup_nog));

		kreg->add(kernel_name::halo_setup, "fortran-msg",
		         boxmg::kernel<grid_topo&,void**>(impls::setup_msg));

		kreg->add(kernel_name::halo_exchange, "fortran-msg",
		          boxmg::kernel<bmg2d::mpi::grid_func&>(impls::msg_exchange));

		kreg->add(kernel_name::halo_stencil_exchange, "fortran-msg",
		          boxmg::kernel<bmg2d::mpi::stencil_op&>(impls::msg_stencil_exchange));

		kreg->add(kernel_name::setup_interp, "fortran-msg",
		          boxmg::kernel<int,int,int,
		          const mpi::stencil_op&,
		          const mpi::stencil_op&,
		          inter::mpi::prolong_op&>(impls::mpi_setup_interp));

		kreg->add(kernel_name::galerkin_prod, "fortran-msg",
		          boxmg::kernel<int,int,int,
		          const inter::mpi::prolong_op&,
		          const mpi::stencil_op&,
		          mpi::stencil_op&>(impls::mpi_galerkin_prod));

		kreg->add(kernel_name::setup_relax, "fortran-msg-rbgs-point",
		         boxmg::kernel<const mpi::stencil_op&,
		         bmg2d::relax_stencil&>(impls::mpi_setup_rbgs_point));

		kreg->add(kernel_name::setup_relax_x, "fortran-msg",
		          boxmg::kernel<const mpi::stencil_op&,
		          bmg2d::relax_stencil&>(impls::mpi_setup_rbgs_x));

		kreg->add(kernel_name::setup_relax_y, "fortran-msg",
		         boxmg::kernel<const mpi::stencil_op&,
		         bmg2d::relax_stencil&>(impls::mpi_setup_rbgs_y));

		kreg->add(kernel_name::setup_cg_boxmg, "fortran-msg",
		          boxmg::kernel<const mpi::stencil_op &,
		          std::shared_ptr<config::reader>,
		          std::shared_ptr<solver>*>(impls::setup_cg_boxmg));

		kreg->add(kernel_name::solve_cg_boxmg, "fortran-msg",
		          boxmg::kernel<const solver &,
		          mpi::grid_func &,
		          const mpi::grid_func &>(impls::solve_cg_boxmg));

		kreg->add(kernel_name::matvec, "fortran-msg",
		          boxmg::kernel<const mpi::stencil_op&,
		          const mpi::grid_func &, mpi::grid_func&>(impls::matvec));

		kreg->add(kernel_name::setup_cg_redist, "c++",
		          boxmg::kernel<const mpi::stencil_op &,
		          std::shared_ptr<config::reader>,
		          std::shared_ptr<mpi::redist_solver>*,
		          std::vector<int>&>(impls::setup_cg_redist));

		kreg->add(kernel_name::solve_cg_redist, "c++",
		          boxmg::kernel<const mpi::redist_solver &,
		          mpi::grid_func &,
		          const mpi::grid_func &>(impls::solve_cg_redist));


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

		return kreg;
	}
}

}}}}
