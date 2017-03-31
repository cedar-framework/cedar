#ifndef CEDAR_PERF_FULL_ITERATOR_H
#define CEDAR_PERF_FULL_ITERATOR_H

#include <cedar/perf/redist_generator.h>


namespace cedar {

class full_iterator : public redist_iterator
{
public:
full_iterator(std::array<int, 2> np, std::array<len_t,2> nglobal,
	                int min_coarse) : redist_iterator(np, nglobal, min_coarse) {
		flag = false;
	}

full_iterator() : redist_iterator() {}

	full_iterator & operator++() {

		if (((nglobal[0] / (nblocks[0]*2)) <= 2*min_coarse) or (np[0] / (nblocks[0]*2) <= 0)) {
			nblocks[1] *= 2;
			nblocks[0] = 1;
		} else {
			nblocks[0] *= 2;
		}

		for (int i = 0; i < 2; i++) {
			nlocal[i] = nglobal[i] / nblocks[i];
		}
		if (!keep_refining(np[0], np[1], nblocks[0], nblocks[1], nlocal[0], nlocal[1], min_coarse)) {
			nblocks[0] = 0; nblocks[1] = 0;
		}

		return *this;
	}

private:
	bool flag;
};


}

#endif
