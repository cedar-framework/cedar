#ifndef CEDAR_2D_MPI_STENCIL_OP_H
#define CEDAR_2D_MPI_STENCIL_OP_H

#include <mpi.h>

#include <cedar/mpi/par_object.h>
#include <cedar/2d/stencil_op.h>

namespace cedar { namespace cdr2 { namespace mpi {

template <class sten>
class stencil_op : public cdr2::stencil_op<sten>, public par_object
{
public:
	stencil_op() {}
stencil_op(topo_ptr grd) :
	cdr2::stencil_op<sten>(grd->nlocal(0)-1, grd->nlocal(1)-1), // remove only one ghost since MPI needs an extra ghost
		par_object(grd, grd->comm)
		{
			// TODO: verify this
			// Fortran kernels expect the extra ghost to be excluded from extents
			this->range_[0] = cedar::range(static_cast<len_t>(1), grd->nlocal(0)-1);
			this->range_[1] = cedar::range(static_cast<len_t>(1), grd->nlocal(1)-1);
			this->len(0)--;
			this->len(1)--;
		}
	using cdr2::stencil_op<sten>::shape;
	using cdr2::stencil_op<sten>::len;
	using cdr2::stencil_op<sten>::set;
	using cdr2::stencil_op<sten>::operator();
	using cdr2::stencil_op<sten>::range;
	using cdr2::stencil_op<sten>::grange;
};

template <class sten>
	std::ostream & operator<<(std::ostream & os, const stencil_op<sten> & obj);

}}}

#endif
