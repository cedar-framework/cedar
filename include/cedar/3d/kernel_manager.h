#ifndef CEDAR_3D_KERNEL_MANAGER_H
#define CEDAR_3D_KERNEL_MANAGER_H

#include <cedar/3d/types.h>
#include <cedar/kernel_manager.h>
#include <cedar/kernel_registry.h>

namespace cedar { namespace cdr3 {

using kernels3 = ser_kernels<stypes>;
using kman_ptr = std::shared_ptr<kernel_manager<kernels3>>;

using point_relax = kernels::point_relax<stypes>;
template<relax_dir rdir>
using plane_relax = kernels::plane_relax<stypes, rdir>;
using coarsen_op = kernels::coarsen_op<stypes>;
using interp_add = kernels::interp_add<stypes>;
using restriction = kernels::restriction<stypes>;
using residual = kernels::residual<stypes>;
using setup_interp = kernels::setup_interp<stypes>;
using solve_cg = kernels::solve_cg<stypes>;

kman_ptr build_kernel_manager(config::reader & conf);
kman_ptr build_kernel_manager(std::shared_ptr<kernel_params> params);

}}

#endif
