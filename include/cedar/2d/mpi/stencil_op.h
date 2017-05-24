#ifndef CEDAR_2D_MPI_STENCIL_OP_H
#define CEDAR_2D_MPI_STENCIL_OP_H

#include <mpi.h>

#include <cedar/mpi/par_object.h>
#include <cedar/2d/stencil_op.h>

namespace cedar { namespace cdr2 { namespace mpi {

template <class sten>
class stencil_op : public stencil_op<sten>, public par_object
{
public:
	stencil_op() {}
stencil_op(topo_ptr grd) :
	stencil_op<sten>(grd->nlocal(0)-1, grd->nlocal(1)-1), // remove only one ghost since MPI needs an extra ghost
		par_object(grd, grd->comm)
		{
			// TODO: verify this
			// Fortran kernels expect the extra ghost to be excluded from extents
			this->len(0)--;
			this->len(1)--;
		}
};

template <class sten>
	std::ostream & operator<<(std::ostream & os, const stencil_op<sten> & obj);

}}}

#endif
