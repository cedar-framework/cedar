#ifndef CEDAR_3D_MPI_STENCIL_OP_H
#define CEDAR_3D_MPI_STENCIL_OP_H

#include <cedar/types.h>
#include <cedar/mpi/par_object.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/mpi/grid_func.h>

namespace cedar { namespace cdr3 { namespace mpi {

template<class sten>
	class stencil_op : public cdr3::stencil_op<sten>, public par_object
{
public:
	stencil_op() {}
stencil_op(topo_ptr grid) :
	cdr3::stencil_op<sten>(grid->nlocal(0)-1, grid->nlocal(1)-1, grid->nlocal(2)-1), // remove only one ghost since MPI needs an extra ghost
		par_object(grid, grid->comm)
		{
			// TODO: verify this
			// Fortran kernels expect the extra ghost to be excluded from extents
			for (auto i : cedar::range<unsigned short>(3)) {
				this->range_[i] = cedar::range(static_cast<len_t>(1), grid->nlocal(i)-1);
				this->grange_[i] = cedar::range(static_cast<len_t>(0), grid->nlocal(i));
				this->len(i)--;
			}
		}

	using cdr3::stencil_op<sten>::shape;
	using cdr3::stencil_op<sten>::len;
	using cdr3::stencil_op<sten>::set;
	using cdr3::stencil_op<sten>::operator();
	using cdr3::stencil_op<sten>::range;
	using cdr3::stencil_op<sten>::grange;
};

template <class sten>
	std::ostream & operator<<(std::ostream & os, const stencil_op<sten> & obj);

}}}

#endif
