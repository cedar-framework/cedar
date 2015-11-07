#include "inter/mpi/prolong_op.h"
#include "restrict_op.h"


namespace boxmg { namespace bmg2d { namespace inter { namespace mpi {

bmg2d::mpi::GridFunc operator*(const RestrictOp & R, const GridFunc &x)
{
	auto & P = R.getP();
	auto & P_par = dynamic_cast<const inter::mpi::ProlongOp&>(P);
	bmg2d::mpi::GridFunc y(P_par.grid_ptr());

	R.apply(x,y);

	return y;
}

}}}}
