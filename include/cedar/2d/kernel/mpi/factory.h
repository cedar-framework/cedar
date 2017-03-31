#ifndef CEDAR_2D_KERNEL_MPI_FACTORY_H
#define CEDAR_2D_KERNEL_MPI_FACTORY_H

#include <memory>

#include <cedar/2d/kernel/mpi/registry.h>

namespace cedar { namespace cdr2 { namespace kernel { namespace mpi {

namespace factory
{
	std::shared_ptr<registry> from_config(config::reader &conf);
	void init(std::shared_ptr<registry> kreg, config::reader & conf);
}

}}}}
#endif
