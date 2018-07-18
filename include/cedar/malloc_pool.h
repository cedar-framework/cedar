#ifndef CEDAR_MALLOC_POOL_H
#define CEDAR_MALLOC_POOL_H

#include <cedar/services/mempool.h>

namespace cedar {

class malloc_pool : public services::mempool
{
public:
	using pool = services::mempool::pool;
	void *addr(memid vtype, std::size_t nbytes) override;
	pool create(memid vtypes, std::size_t nbytes) override;
	int pos(std::size_t nbytes) override;

	~malloc_pool();

protected:
	std::vector<void*> addrs;
};

}
#endif
