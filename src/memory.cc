#include <cstdlib>
#include <iostream>

#ifdef CUDA_MANAGED
#include <cuda_runtime.h>
#endif

#include <cedar/memory.h>

namespace cedar { namespace memory {

std::function<void*(std::size_t nbytes)> allocator = std::malloc;
std::function<void(void * addr)> deallocator = std::free;

void init(config & conf)
{
	auto mtype = conf.get<std::string>("memory", "malloc");
	#ifdef CUDA_MANAGED
	if (mtype == "managed") {
		allocator = [](std::size_t nbytes)
		            {
			            void *addr;
			            cudaMallocManaged(&addr, nbytes, cudaMemAttachGlobal);
			            return addr;
		            };
		deallocator = cudaFree;
	}
	#endif
}

}}
