#include <cstdlib>
#include <cedar/malloc_pool.h>

namespace cedar {

	void * malloc_pool::addr(vecs vtype, std::size_t nbytes)
	{
		void *ret = std::malloc(nbytes);

		addrs.push_back(ret);
		return ret;
	}


	malloc_pool::~malloc_pool()
	{
		for (std::size_t i = 0; i < addrs.size(); i++) {
			std::free(addrs[i]);
		}
	}

}

