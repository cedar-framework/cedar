#ifndef CEDAR_2D_INTER_MPI_PROLONG_OP_H
#define CEDAR_2D_INTER_MPI_PROLONG_OP_H

#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/inter/types.h>


namespace cedar { namespace cdr2 { namespace inter { namespace mpi {


namespace mpi = cedar::cdr2::mpi;

class prolong_op : public cdr2::stencil_op<inter::dir>, public par_object
{
public:
	prolong_op() {};
	prolong_op(topo_ptr grid);
	friend std::ostream & operator<< (std::ostream &os, const prolong_op & P);

	using cdr2::stencil_op<inter::dir>::shape;
	using cdr2::stencil_op<inter::dir>::len;
	using cdr2::stencil_op<inter::dir>::set;
	using cdr2::stencil_op<inter::dir>::operator();
	using cdr2::stencil_op<inter::dir>::range;
	using cdr2::stencil_op<inter::dir>::grange;

	mpi::stencil_op<five_pt> * fine_op_five;
	mpi::stencil_op<nine_pt> * fine_op_nine;
	bool fine_is_five;
	mpi::grid_func *residual;
};

}}}}


#endif
