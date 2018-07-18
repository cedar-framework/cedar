#include <cedar/services/mempool.h>

extern "C" {
	void cedar_mempool_pos(cedar::services::mempool *mp,
	                       int nbytes, int *pos)
	{
		std::size_t nbytes_in;
		nbytes_in = nbytes;

		*pos = mp->pos(nbytes_in);
	}
}
