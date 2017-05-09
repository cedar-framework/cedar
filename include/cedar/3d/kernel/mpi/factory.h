#ifndef CEDAR_3D_KERNEL_MPI_FACTORY_H
#define CEDAR_3D_KERNEL_MPI_FACTORY_H

#include <memory>

#include <cedar/3d/kernel/mpi/registry.h>

namespace cedar { namespace cdr3 { namespace kernel { namespace mpi {

namespace factory
{
	std::shared_ptr<registry> from_config(config::reader &conf);
	void init(std::shared_ptr<registry> kreg, config::reader & conf);
}

}}}}

#endif
