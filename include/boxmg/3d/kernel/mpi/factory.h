#ifndef BOXMG_3D_KERNEL_MPI_FACTORY_H
#define BOXMG_3D_KERNEL_MPI_FACTORY_H

#include <memory>

#include <boxmg/3d/kernel/mpi/registry.h>

namespace boxmg { namespace bmg3 { namespace kernel { namespace mpi {

namespace factory
{
	std::shared_ptr<registry> from_config(config::Reader &conf);
}

}}}}

#endif
