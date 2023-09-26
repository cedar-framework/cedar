#ifndef CEDAR_DEVICE_H
#define CEDAR_DEVICE_H

#include <ftl/Base.hpp>

namespace cedar {
    using cpu = ftl::device::CPU;
    using gpu = ftl::device::GPU;
    using alloc_device = ftl::BufferAllocateDevice;
}

#endif
