#ifndef BOXMG_PARTITION_H
#define BOXMG_PARTITION_H

#include <boxmg/types.h>

namespace boxmg {

class partition
{
public:
partition(len_t n, len_t num_procs) : num_procs(num_procs), n(n) {}
	virtual ~partition() {};
	virtual len_t low(len_t id) = 0;
	virtual len_t high(len_t id) = 0;
	virtual len_t size(len_t id) = 0;
	virtual len_t owner(len_t index) = 0;

protected:
	len_t num_procs;
	len_t n;
};

}
#endif
