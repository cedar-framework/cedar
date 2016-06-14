#ifndef BOXMG_PERF_REDIST_GENERATOR_H
#define BOXMG_PERF_REDIST_GENERATOR_H

#include <boxmg/types.h>
#include <array>

namespace boxmg {

static inline bool keep_refining(int npx, int npy, int nbx, int nby, len_t nlx, len_t nly,
                          int min_coarse)
{
	bool ret = ((npx / nbx) > 0 and (npy / nby) > 0);
	//ret = ret and ((npx / nbx) * (npy / nby) > 1);
	ret = ret and (nlx > 2*min_coarse);
	ret = ret and (nly > 2*min_coarse);

	return ret;
}


class redist_iterator
{
public:
	typedef std::array<int,2> value_type;

redist_iterator(std::array<int, 2> np, std::array<len_t,2> nglobal,
                int min_coarse):
	nblocks({1,1}), nglobal(nglobal), nlocal(nglobal), np(np), min_coarse(min_coarse) {}
redist_iterator() : nblocks({0,0}) {}

	value_type operator*() {
		return nblocks;
	}

	redist_iterator & operator++() {
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

	friend bool operator==(const redist_iterator & lhs, const redist_iterator & rhs)
	{
		return lhs.nblocks == rhs.nblocks;
	}

	friend bool operator!=(const redist_iterator & lhs, const redist_iterator & rhs)
	{
		return ! (lhs == rhs);
	}


private:
	std::array<int,2> nblocks;
	std::array<len_t,2> nglobal;
	std::array<len_t,2> nlocal;
	std::array<int,2> np;
	int min_coarse;
};


class redist_generator
{
public:
	redist_generator(std::array<int, 2> np, std::array<len_t, 2> nglobal,
	                 int min_coarse) :
	np(np), nglobal(nglobal), min_coarse(min_coarse) {}

	redist_iterator begin() {
		return redist_iterator(np, nglobal, min_coarse);
	}


	redist_iterator end() {
		return redist_iterator();
	}

private:
	std::array<int,2> np;
	std::array<len_t,2> nglobal;
	int min_coarse;
};

}
#endif
