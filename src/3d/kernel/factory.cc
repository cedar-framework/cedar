#include <cedar/kernel.h>
#include <cedar/kernel_name.h>
#include <cedar/3d/residual.h>
#include <cedar/3d/relax/setup_relax.h>
#include <cedar/3d/inter/setup_interp.h>
#include <cedar/3d/inter/galerkin_prod.h>
#include <cedar/3d/relax/relax.h>
#include <cedar/3d/cg/setup_cg_lu.h>
#include <cedar/3d/cg/solve_cg.h>
#include <cedar/3d/inter/restrict.h>
#include <cedar/3d/inter/interp.h>

#include <cedar/3d/kernel/factory.h>

namespace cedar { namespace cdr3 { namespace kernel {

namespace factory
{
	namespace name = cedar::kernel_name;

	std::shared_ptr<registry> from_config(config::reader & conf)
	{
		auto kreg = std::make_shared<registry>();

		auto params = build_kernel_params(conf);

		kreg->add(name::residual, "fortran",
		          cedar::kernel<const stencil_op &,
		          const grid_func &,
		          const grid_func &,
		          grid_func&>(impls::residual, params));

		kreg->add(name::setup_interp, "fortran",
		          cedar::kernel<int, int, int,
		          const stencil_op &,
		          const stencil_op &,
		          inter::prolong_op &>(impls::setup_interp, params));

		kreg->add(name::setup_relax, "fortran-rbgs-point",
		          cedar::kernel<const stencil_op&,
		          relax_stencil&>(impls::setup_rbgs_point, params));

		// kreg->add(name::setup_relax_xy, "c++",
		//           cedar::kernel<const stencil_op &,
		//           relax_stencil&>(impls::setup_relax_xy));

		kreg->add(name::galerkin_prod, "fortran-ex",
		         cedar::kernel<int,int,int,
		         const inter::prolong_op&,
		         const stencil_op&,
		          stencil_op&>(impls::galerkin_prod, params));

		kreg->add(name::relax, "fortran-rbgs-point",
		          cedar::kernel<const stencil_op&,
		          grid_func&,
		          const grid_func&,
		          const relax_stencil&,
		          cycle::Dir>(impls::relax_rbgs_point, params));

		kreg->add(name::setup_cg_lu, "fortran",
		         cedar::kernel<const stencil_op&,
		          grid_func&>(impls::setup_cg_lu, params));

		kreg->add(name::solve_cg, "fortran",
		          cedar::kernel<grid_func&,
		          const grid_func&,
		          const grid_func&,
		          real_t*>(impls::fortran_solve_cg, params));

		kreg->add(name::restriction, "fortran",
		         cedar::kernel<const inter::restrict_op&,
		         const grid_func&,
		          grid_func&>(impls::fortran_restrict, params));

		kreg->add(name::interp_add,"fortran",
		         cedar::kernel<const inter::prolong_op&,
		         const grid_func&,
		         const grid_func&,
		          grid_func&>(impls::fortran_interp, params));

		std::vector<std::tuple<std::string, std::string, std::string>> defaults = {
			std::make_tuple(name::residual, "kernels.residual", "fortran"),
			std::make_tuple(name::galerkin_prod, "kernels.galerkin-prod", "fortran-ex"),
			std::make_tuple(name::setup_interp, "kernels.setup-interp", "fortran"),
			std::make_tuple(name::setup_relax, "kernels.setup-relax", "fortran-rbgs-point"),
			// std::make_tuple(name::setup_relax_xy, "kernels.setup-relax-xy", "c++"),
			std::make_tuple(name::restriction, "kernels.restrict", "fortran"),
			std::make_tuple(name::relax, "kernels.relax", "fortran-rbgs-point"),
			std::make_tuple(name::setup_cg_lu, "kernels.setup-cg-lu", "fortran"),
			std::make_tuple(name::interp_add, "kernels.interp-add", "fortran"),
			std::make_tuple(name::solve_cg, "kernels.solve-cg", "fortran")
		};

		for (auto&& v : defaults) {
			std::string kname = conf.get<std::string>(std::get<1>(v), std::get<2>(v));
			log::info << "Using '" + kname + "' for " <<  std::get<0>(v) << "." << std::endl;
			kreg->set(std::get<0>(v), kname);
		}

		return kreg;
	}
}

}}}
