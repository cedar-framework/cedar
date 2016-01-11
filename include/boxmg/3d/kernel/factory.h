#ifndef BOXMG_3D_KERNEL_FACTORY_H
#define BOXMG_3D_KERNEL_FACTORY_H

#include <memory>

#include <boxmg/3d/kernel/registry.h>

namespace boxmg { namespace bmg3 { namespace kernel {

namespace factory
{
	std::shared_ptr<registry> from_config(config::Reader &conf);
}

}}}

#endif
