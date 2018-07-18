#ifndef CEDAR_3D_MPI_MEMPOOL_H
#define CEDAR_3D_MPI_MEMPOOL_H

#include <array>
#include <map>

#include <cedar/services/mempool.h>

namespace cedar { namespace cdr3 { namespace mpi {

class plane_mempool : public services::mempool
{
public:
	using pool = services::mempool::pool;
	using addrs_type = std::array<std::map<std::size_t, char*>, num_memid>;
	plane_mempool(int nplanes);
	plane_mempool(int worker_id, addrs_type *addrs);
	addrs_type *get_addrs() { return addrs; }
	void *addr(memid vtype, std::size_t nbytes) override;
	pool create(memid vtype, std::size_t nbytes) override;
	int pos(std::size_t nbytes) override;
	~plane_mempool();

protected:
	int wid;                  /** worker id (index of plane) */
	bool ismaster;
	int nplanes;
	addrs_type *addrs;
	std::unique_ptr<addrs_type> addrs_data;
};


}}}

#endif
