#ifndef CEDAR_MEMPOOL_H
#define CEDAR_MEMPOOL_H

#include <cedar/service.h>

namespace cedar { namespace services {

class mempool : public service
{
public:
	enum vecs { rhs, sol, res };
	const static std::string name() { return "memory pool"; }
	virtual void *addr(enum vecs vtype, std::size_t nbytes) = 0;
};

}}

#endif
