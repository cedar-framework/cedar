#ifndef CEDAR_2D_KERNEL_FACTORY_H
#define CEDAR_2D_KERNEL_FACTORY_H

#include <memory>

#include <cedar/2d/kernel/registry.h>

namespace cedar { namespace cdr2 { namespace kernel {

namespace factory
{
	std::shared_ptr<registry> from_config(config::reader &conf);
}

}}}
#endif
