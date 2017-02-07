#ifndef BOXMG_PERF_REDIST_GENERATOR_H
#define BOXMG_PERF_REDIST_GENERATOR_H

#include <boxmg/types.h>
#include <array>

namespace boxmg {

static inline bool keep_refining(int npx, int npy, int nbx, int nby, len_t nlx, len_t nly,
                                 len_t min_coarse)
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
                len_t min_coarse):
	nblocks({1,1}), nglobal(nglobal), nlocal(nglobal), np(np), min_coarse(min_coarse) {}
redist_iterator() : nblocks({0,0}) {}

	value_type operator*() {
		return nblocks;
	}

	virtual redist_iterator & operator++() = 0;

	friend bool operator==(const redist_iterator & lhs, const redist_iterator & rhs)
	{
		return lhs.nblocks == rhs.nblocks;
	}

	friend bool operator!=(const redist_iterator & lhs, const redist_iterator & rhs)
	{
		return ! (lhs == rhs);
	}


protected:
	std::array<int,2> nblocks;
	std::array<len_t,2> nglobal;
	std::array<len_t,2> nlocal;
	std::array<int,2> np;
	len_t min_coarse;
};


template <class rit>
class redist_generator
{
public:
	redist_generator(std::array<int, 2> np, std::array<len_t, 2> nglobal,
	                 len_t min_coarse) :
	np(np), nglobal(nglobal), min_coarse(min_coarse) {}

	rit begin() {
		return rit(np, nglobal, min_coarse);
	}


	rit end() {
		return rit();
	}

private:
	std::array<int,2> np;
	std::array<len_t,2> nglobal;
	len_t min_coarse;
};

}
#endif
