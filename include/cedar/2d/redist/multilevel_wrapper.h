#ifndef CEDAR_2D_REDIST_MULTILEVEL_WRAPPER_H
#define CEDAR_2D_REDIST_MULTILEVEL_WRAPPER_H

#include <cedar/2d/types.h>
#include <cedar/2d/mpi/types.h>
#include <cedar/config.h>
#include <cedar/2d/mpi/msg_exchanger.h>

namespace cedar { namespace cdr2 {

template<class solver>
struct serial_type { static const bool value; };


		template<class inner_solver>
		class multilevel_wrapper
		{
		public:
			static const bool is_serial = serial_type<inner_solver>::value;
			using stypes = typename std::conditional<is_serial, cdr2::stypes, mpi::stypes>::type;
			using stencil_op = typename stypes::template stencil_op<nine_pt>;
			using grid_func = typename stypes::grid_func;
		multilevel_wrapper(stencil_op & sop) : inner(sop) {}
		multilevel_wrapper(stencil_op & sop, std::shared_ptr<config> conf) : inner(sop, conf) {}
			void cycle(grid_func & x, const grid_func & b)
			{
				inner.vcycle(x, b);
			}

			config & get_config() { return inner.get_config(); }

		protected:
			inner_solver inner;
		};

		template<class fsten> class solver;
		namespace mpi {
			template<class fsten> class solver;
		}
template<>
struct serial_type<cdr2::solver<nine_pt>> { static const bool value = true; };

template<>
struct serial_type<cdr2::mpi::solver<nine_pt>> { static const bool value = false; };
}}
#endif
