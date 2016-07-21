#ifndef BOXMG_PERF_GREEDY_ITERATOR_H
#define BOXMG_PERF_GREEDY_ITERATOR_H

#include <boxmg/perf/redist_generator.h>


namespace boxmg {

class greedy_iterator : public redist_iterator
{
public:
greedy_iterator(std::array<int, 2> np, std::array<len_t,2> nglobal,
	                int min_coarse) : redist_iterator(np, nglobal, min_coarse) {}

greedy_iterator() : redist_iterator() {}

	greedy_iterator & operator++() {
		if (nlocal[0] > nlocal[1]) {
			if (((nglobal[0] / (nblocks[0]*2)) <= 2*min_coarse) or (np[0] / (nblocks[0]*2) <= 0))
				nblocks[1] *= 2;
			else
				nblocks[0] *= 2;
		} else {
			if (((nglobal[1] / (nblocks[1]*2)) <= 2*min_coarse) or (np[1] / (nblocks[1]*2) <=0))
				nblocks[0] *= 2;
			else
				nblocks[1] *= 2;
		}
		for (int i = 0; i < 2; i++) {
			nlocal[i] = nglobal[i] / nblocks[i];
		}
		if (!keep_refining(np[0], np[1], nblocks[0], nblocks[1], nlocal[0], nlocal[1], min_coarse)) {
			nblocks[0] = 0; nblocks[1] = 0;
		}
	}
};


}

#endif
