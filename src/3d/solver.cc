#include <algorithm>

#include <boxmg/3d/kernel/factory.h>
#include <boxmg/3d/kernel/registry.h>
#include <boxmg/3d/solver.h>

using namespace boxmg;
using namespace boxmg::bmg3;

solver::solver(stencil_op&& fop)
{


}


void solver::add_level(stencil_op & fop, int num_levels)
{

}


int solver::compute_num_levels(stencil_op & fop)
{


}


std::shared_ptr<bmg3::kernel::registry> solver::kernel_registry()
{
	return std::static_pointer_cast<bmg3::kernel::registry>(kreg);
}
