#ifndef CEDAR_HALO_EXCHANGER_H
#define CEDAR_HALO_EXCHANGER_H

#include <cedar/types.h>

namespace cedar {

	class halo_exchanger
	{
	public:
		virtual void exchange_func(int k, real_t *gf) = 0;
		virtual void exchange_sten(int k, real_t *so) = 0;
	};
}

#endif
