#ifndef BOXMG_2D_KERNEL_FACTORY_H
#define BOXMG_2D_KERNEL_FACTORY_H

#include <memory>

#include <boxmg/2d/kernel/registry.h>

namespace boxmg { namespace bmg2d { namespace kernel {

namespace factory
{
	std::shared_ptr<registry> from_config(config::reader &conf);
}

}}}
#endif
