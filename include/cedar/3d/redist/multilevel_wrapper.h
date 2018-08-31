#ifndef CEDAR_3D_REDIST_MULTILEVEL_WRAPPER_H
#define CEDAR_3D_REDIST_MULTILEVEL_WRAPPER_H

#include <cedar/config.h>
#include <cedar/3d/types.h>
#include <cedar/3d/mpi/types.h>
#include <cedar/3d/mpi/msg_exchanger.h>

namespace cedar { namespace cdr3 {

template<class solver>
struct serial_type { static const bool value; };


		template<class inner_solver>
		class multilevel_wrapper
		{
		public:
			static const bool is_serial = serial_type<inner_solver>::value;
			using stypes = typename std::conditional<is_serial, cdr3::stypes, mpi::stypes>::type;
			using stencil_op = typename stypes::template stencil_op<xxvii_pt>;
			using grid_func = typename stypes::grid_func;
		multilevel_wrapper(stencil_op & sop) : inner(sop) {}
		multilevel_wrapper(stencil_op & sop, std::shared_ptr<config> conf) : inner(sop, conf) {}
			void cycle(grid_func & x, const grid_func & b)
			{
				inner.vcycle(x, b);
			}

			config & get_config() { return inner.get_config(); }

			service_manager<mpi::stypes> & get_services() { return inner.get_kernels()->services(); }

			inner_solver & get_inner() { return inner; }

		protected:
			inner_solver inner;
		};

		template<class fsten> class solver;
		namespace mpi {
		template<class sten>
			class solver;
		}

template<>
struct serial_type<cdr3::solver<xxvii_pt>> { static const bool value = true; };

template<>
struct serial_type<cdr3::mpi::solver<xxvii_pt>> { static const bool value = false; };

}}

#endif
