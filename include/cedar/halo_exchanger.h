#ifndef CEDAR_HALO_EXCHANGER_H
#define CEDAR_HALO_EXCHANGER_H

#include <cedar/types.h>

namespace cedar {
	template<short ND>
	struct halo_exchanger
	{
		std::function<void(int k, int nog, real_t*, std::array<len_t, ND> len, void*)> exchange;
		std::function<void(int,int,real_t*,void*)> stencil_exchange;
	};
}

#endif
