#include <memory>
#include <boxmg/kernel.h>
#include <boxmg/kernel_name.h>
#include <boxmg/3d/mpi/halo.h>
#include <boxmg/3d/residual.h>

#include <boxmg/3d/kernel/mpi/factory.h>

namespace boxmg { namespace bmg3 { namespace kernel { namespace mpi {

namespace mpi = boxmg::bmg3::mpi;

namespace factory
{
	std::shared_ptr<registry> from_config(config::Reader &conf)
	{
		auto kreg = std::make_shared<registry>();

		kreg->add(kernel_name::halo_setup, "fortran-msg",
		         boxmg::kernel<grid_topo&,void**>(impls::setup_msg));

		kreg->add(kernel_name::halo_exchange, "fortran-msg",
		          boxmg::kernel<bmg3::mpi::grid_func&>(impls::msg_exchange));

		kreg->add(kernel_name::halo_stencil_exchange, "fortran-msg",
		          boxmg::kernel<bmg3::mpi::stencil_op&>(impls::msg_stencil_exchange));

		kreg->add(kernel_name::residual, "fortran-msg",
		         boxmg::kernel<const mpi::stencil_op &,
		          const mpi::grid_func &,
		          const mpi::grid_func &,
		          mpi::grid_func&>(impls::mpi_residual_fortran));

		std::vector<std::tuple<std::string, std::string, std::string>> defaults = {
			std::make_tuple(kernel_name::residual, "kernels.residual", "fortran-msg"),
			std::make_tuple(kernel_name::halo_setup, "kernels.halo-setup", "fortran-msg"),
			std::make_tuple(kernel_name::halo_exchange, "kernels.halo-exchange", "fortran-msg"),
			std::make_tuple(kernel_name::halo_stencil_exchange, "kernels.halo-stencil-exchange", "fortran-msg")
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
