#ifndef CEDAR_3D_INTER_MPI_PROLONG_OP_H
#define CEDAR_3D_INTER_MPI_PROLONG_OP_H

#include <cedar/3d/inter/types.h>
#include <cedar/3d/stencil_op.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/grid_func.h>

namespace cedar { namespace cdr3 { namespace inter { namespace mpi {

namespace mpi = cedar::cdr3::mpi;

class prolong_op : public cdr3::stencil_op<inter::dir>, public par_object
{
public:
	prolong_op() {};
	prolong_op(topo_ptr grid);
	friend std::ostream & operator<< (std::ostream &os, const prolong_op & P);

	using cdr3::stencil_op<inter::dir>::shape;
	using cdr3::stencil_op<inter::dir>::len;
	using cdr3::stencil_op<inter::dir>::set;
	using cdr3::stencil_op<inter::dir>::operator();
	using cdr3::stencil_op<inter::dir>::range;
	using cdr3::stencil_op<inter::dir>::grange;

	mpi::stencil_op<seven_pt> * fine_op_seven;
	mpi::stencil_op<xxvii_pt> * fine_op_xxvii;
	bool fine_is_seven;
};

}}}}

#endif
