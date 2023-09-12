#ifndef CEDAR_2D_GPU_KERNEL_MANAGER_H
#define CEDAR_2D_GPU_KERNEL_MANAGER_H

#include <cedar/2d/gpu/types.h>
#include <cedar/kernel_manager.h>
#include <cedar/kernel_registry.h>
#include <cedar/service_manager.h>

namespace cedar { namespace cdr2 { namespace gpu { namespace mpi {

using kman_ptr = std::shared_ptr<kernel_manager<klist<stypes, exec_mode::mpi>, stypes>>;

using point_relax = kernels::point_relax<stypes>;
using coarsen_op = kernels::coarsen_op<stypes>;
using interp_add = kernels::interp_add<stypes>;
using restriction = kernels::restriction<stypes>;
using matvec = kernels::matvec<stypes>;
using residual = kernels::residual<stypes>;
using setup_interp = kernels::setup_interp<stypes>;
using setup_nog = kernels::setup_nog<stypes>;
using solve_cg = kernels::solve_cg<stypes>;
using halo_exchange = services::halo_exchange<stypes>;
using message_passing = services::message_passing;
using mempool = services::mempool;

kman_ptr build_kernel_manager(config & conf);
kman_ptr build_kernel_manager(std::shared_ptr<kernel_params> params);

}}}}
#endif
