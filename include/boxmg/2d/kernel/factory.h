#ifndef BOXMG_2D_KERNEL_FACTORY_H
#define BOXMG_2D_KERNEL_FACTORY_H

#include <memory>

#include "registry.h"

namespace boxmg { namespace bmg2d { namespace kernel {

namespace factory
{
	std::shared_ptr<Registry> from_config(config::Reader &conf);
}

}}}
#endif
