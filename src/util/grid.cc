#include <cedar/util/grid.h>

namespace cedar
{

	AlignedVector<real_t> linspace(real_t start, real_t end, len_t num)
	{
		AlignedVector<real_t> v(num);

		real_t step = (end-start)/(num-1);

		for (auto i: range(num)) v[i] = start + i*step;

		return v;
	}

}
