#ifndef BOXMG_2D_KERNEL_MPI_FACTORY_H
#define BOXMG_2D_KERNEL_MPI_FACTORY_H

#include <memory>

#include <boxmg/2d/kernel/mpi/registry.h>

namespace boxmg { namespace bmg2d { namespace kernel { namespace mpi {

namespace factory
{
	std::shared_ptr<registry> from_config(config::reader &conf);
	void init(std::shared_ptr<registry> kreg, config::reader & conf);
}

}}}}
#endif
