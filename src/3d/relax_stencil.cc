#include <vector>

#include <cedar/3d/relax_stencil.h>

using namespace cedar;
using namespace cedar::cdr3;

relax_stencil::relax_stencil(len_t nx, len_t ny, len_t nz, unsigned int nghosts):
	array<real_t,4>(nx+3*nghosts, ny + 3*nghosts, nz + 3*nghosts, 2)
{
	num_ghosts = nghosts;
	std::vector<len_t> lens = {nx, ny, nz};

	for (auto i : cedar::range(3)) {
		range_[i] = cedar::range(static_cast<len_t>(nghosts), static_cast<len_t>(lens[i] + nghosts));
		grange_[i] = cedar::range(static_cast<len_t>(0), lens[i] + 3*nghosts);
	}
}
