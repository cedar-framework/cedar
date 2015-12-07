#include <boxmg/2d/inter/mpi/prolong_op.h>
#include <boxmg/2d/inter/restrict_op.h>
#include <boxmg/2d/mpi/grid_func.h>


namespace boxmg { namespace bmg2d { namespace inter { namespace mpi {

bmg2d::mpi::grid_func operator*(const restrict_op & R, const grid_func &x)
{
	auto & P = R.getP();
	auto & P_par = dynamic_cast<const inter::mpi::prolong_op&>(P);
	bmg2d::mpi::grid_func y(P_par.grid_ptr());

	R.apply(x,y);

	return y;
}

}}}}
