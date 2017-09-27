#ifndef CEDAR_HALO_EXCHANGER_H
#define CEDAR_HALO_EXCHANGER_H

#include <cedar/types.h>
#include <cedar/3d/mpi/msg_exchanger.h>
#include <cedar/2d/mpi/msg_exchanger.h>


namespace cedar {

	template <unsigned int nd>
		using halo_exchanger = typename std::conditional<nd==2,
		cdr2::mpi::msg_exchanger,
		cdr3::mpi::msg_exchanger>::type;
}

#endif
