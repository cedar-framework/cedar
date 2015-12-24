#ifndef BOXMG_KERNEL_NAME
#define BOXMG_KERNEL_NAME

#include <string>
#include <tuple>


namespace boxmg {

namespace kernel_name {
	const std::string residual("residual");
	const std::string setup_interp("setup-oi-interp");
	const std::string galerkin_prod("galerkin-prod");
	const std::string setup_relax("setup pointwise relaxation");
	const std::string setup_relax_x("setup x-line relaxation");
	const std::string setup_relax_y("setup y-line relaxation");
	const std::string setup_cg_lu("setup-cg-lu");
	const std::string relax("pointwise relaxation");
	const std::string relax_lines_x("x-line relaxation");
	const std::string relax_lines_y("y-line relaxation");
	const std::string restriction("restriction");
	const std::string interp_add("interpolate and add");
	const std::string solve_cg("coarse grid solve");
	const std::string setup_nog("setup number of grids");
	const std::string halo_setup("halo exchange setup");
	const std::string halo_exchange("halo exchange");
	const std::string halo_stencil_exchange("halo stencil exchange");
	const std::string setup_cg_boxmg("BoxMG coarse grid solve setup");
	const std::string solve_cg_boxmg("BoxMG coarse grid solver");
	const std::string matvec("matrix vector multiply");
}

}
#endif
