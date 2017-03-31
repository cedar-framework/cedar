#ifndef CEDAR_BLOCKPARTITION_H
#define CEDAR_BLOCKPARTITION_H

#include <cedar/mpi/partition.h>

namespace cedar
{
	class block_partition: public partition
	{
	public:
	block_partition(len_t n, len_t num_procs):
		partition(n, num_procs) {};

		inline len_t low(len_t id)
		{
			return id * n / num_procs;
		};

		inline len_t high(len_t id)
		{
			return low(id+1) - 1;
		};

		inline len_t size(len_t id)
		{
			return high(id) - low(id) + 1;
		};

		inline len_t owner(len_t index)
		{
			return (num_procs*(index+1)-1)/n;
		};
	};
}

#endif
