#include <cstdlib>
#include <cedar/malloc_pool.h>

namespace cedar {

	void * malloc_pool::addr(memid vtype, std::size_t nbytes)
	{
		void *ret = std::malloc(nbytes);

		addrs.push_back(ret);
		return ret;
	}


	malloc_pool::pool malloc_pool::create(memid vtype, std::size_t nbytes)
	{
		pool ret;

		ret.addr = (char*) std::malloc(nbytes);
		ret.size = nbytes;
		ret.fullsize = nbytes;

		addrs.push_back((void*) ret.addr);

		return ret;
	}


	int malloc_pool::pos(std::size_t nbytes) { return 0; }


	malloc_pool::~malloc_pool()
	{
		for (std::size_t i = 0; i < addrs.size(); i++) {
			std::free(addrs[i]);
		}
	}

}

