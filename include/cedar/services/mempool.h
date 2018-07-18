#ifndef CEDAR_MEMPOOL_H
#define CEDAR_MEMPOOL_H

#include <cedar/service.h>

namespace cedar { namespace services {

class mempool : public service
{
public:
	struct pool
	{
		int fullsize;
		int size;
		char *addr;
	};
	enum memid { rhs, sol, res, x_redist, b_redist,
	             tricomm_group_x, tricomm_group_y, tricomm_iface_x, tricomm_iface_y, num_memid };
	const static std::string name() { return "memory pool"; }
	virtual void *addr(memid vtype, std::size_t nbytes) = 0;
	virtual pool create(memid vtype, std::size_t nbytes) = 0;
	virtual int pos(std::size_t nbytes) = 0;
};

}}

#endif
