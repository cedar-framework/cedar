#include <cstdlib>

#include <cedar/types.h>
#include <cedar/3d/mpi/plane_mempool.h>


using namespace cedar;
using namespace cedar::cdr3;
using namespace cedar::cdr3::mpi;

plane_mempool::plane_mempool(int nplanes) : wid(0), ismaster(true), nplanes(nplanes)
{
	addrs_data = std::make_unique<addrs_type>();
	addrs = addrs_data.get();
}


plane_mempool::plane_mempool(int worker_id, addrs_type *addrs) :
	wid(worker_id), ismaster(false), nplanes(-1), addrs(addrs) {}


void *plane_mempool::addr(memid vtype, std::size_t nbytes)
{
	char *mloc;

	auto & amap = (*addrs)[vtype];
	if (ismaster) {
		amap[nbytes] = (char*) std::malloc(nbytes*nplanes);
		mloc = amap[nbytes];
	} else {
		mloc = amap[nbytes];
		mloc += wid * nbytes;
	}

	return mloc;
}


plane_mempool::pool plane_mempool::create(memid vtype, std::size_t nbytes)
{
	pool ret;

	auto &amap = (*addrs)[vtype];
	if (ismaster)
		amap[nbytes] = (char*) std::malloc(nbytes*nplanes);
	ret.addr = amap[nbytes];
	ret.fullsize = nbytes*nplanes;
	ret.size = nbytes;

	return ret;
}


int plane_mempool::pos(std::size_t nbytes)
{
	return wid * nbytes;
}


plane_mempool::~plane_mempool()
{
	if (ismaster) {
		for (auto & amap : (*addrs)) {
			for (auto &kv : amap) {
				free(kv.second);
			}
		}
	}
}
