#ifndef CEDAR_HALO_EXCHANGER_H
#define CEDAR_HALO_EXCHANGER_H

#include <cedar/types.h>
#include <cedar/3d/mpi/msg_exchanger.h>
#include <cedar/2d/mpi/msg_exchanger.h>


namespace cedar {

	template<unsigned int nd>
		struct halo_exchanger_t {
			typedef void type;
		};

	template<>
		struct halo_exchanger_t<2> {
		typedef cedar::cdr2::mpi::msg_exchanger type;
	};

	template<>
		struct halo_exchanger_t<3> {
		typedef cdr3::mpi::msg_exchanger type;
	};

	template <unsigned int nd>
		using halo_exchanger = typename halo_exchanger_t<nd>::type;
}

#endif
