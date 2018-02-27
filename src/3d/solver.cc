#include <algorithm>

// TODO: remove include
// #include <cedar/3d/relax/setup_planes.h>
// #include <cedar/3d/relax/relax_planes.h>

#include <cedar/3d/kernel/registry.h>
#include <cedar/3d/solver.h>

using namespace cedar;
using namespace cedar::cdr3;

// solver::~solver()
// {
// 	delete[] bbd;
// }

// void solver::setup_relax_plane(stencil_op & sop, bmg_level & level)
// {
// 	kernel::impls::setup_relax_xy(sop, level.planes_xy);
// 	kernel::impls::setup_relax_xz(sop, level.planes_xz);
// 	kernel::impls::setup_relax_yz(sop, level.planes_yz);
// }


// void solver::relax_plane(const stencil_op & so, grid_func & x,
//                          const grid_func & b, cycle::Dir cdir,
//                          bmg_level & level)
// {
// 	if (cdir == cycle::Dir::DOWN) {
// 		kernel::impls::relax_xy(so, x, b, cdir, level.planes_xy);
// 		kernel::impls::relax_yz(so, x, b, cdir, level.planes_yz);
// 		kernel::impls::relax_xz(so, x, b, cdir, level.planes_xz);
// 	} else {
// 		kernel::impls::relax_xz(so, x, b, cdir, level.planes_xz);
// 		kernel::impls::relax_yz(so, x, b, cdir, level.planes_yz);
// 		kernel::impls::relax_xy(so, x, b, cdir, level.planes_xy);
// 	}

// }
