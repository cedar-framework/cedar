#ifndef CEDAR_3D_INTER_MPI_PROLONG_OP_H
#define CEDAR_3D_INTER_MPI_PROLONG_OP_H

#include <cedar/3d/stencil_op.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/mpi/grid_func.h>
#include <cedar/3d/prolong_op.h>

namespace cedar { namespace cdr3 { namespace mpi {


class prolong_op : public cdr3::stencil_op<inter_dir>, public par_object
{
public:
	prolong_op() {};
	prolong_op(topo_ptr grid);
	friend std::ostream & operator<< (std::ostream &os, const prolong_op & P);

	using cdr3::stencil_op<inter_dir>::shape;
	using cdr3::stencil_op<inter_dir>::len;
	using cdr3::stencil_op<inter_dir>::set;
	using cdr3::stencil_op<inter_dir>::operator();
	using cdr3::stencil_op<inter_dir>::range;
	using cdr3::stencil_op<inter_dir>::grange;

	mpi::stencil_op<seven_pt> * fine_op_seven;
	mpi::stencil_op<xxvii_pt> * fine_op_xxvii;
	bool fine_is_seven;
};

}}}

#endif
