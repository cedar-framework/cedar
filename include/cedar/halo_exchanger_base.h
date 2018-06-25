#ifndef CEDAR_HALO_EXCHANGER_BASE_H
#define CEDAR_HALO_EXCHANGER_BASE_H

#include <cedar/array.h>
#include <cedar/types.h>

namespace cedar {

	class halo_exchanger_base
	{
	public:
		virtual void exchange_func(int k, real_t *gf) = 0;
		virtual void exchange_sten(int k, real_t *so) = 0;
		virtual aarray<int, len_t, 2> & leveldims(int k) = 0;
		virtual len_t * datadist(int k, int grid) = 0;
		virtual std::vector<real_t> & linebuf() = 0;
	};
}

#endif
