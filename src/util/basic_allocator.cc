#include <cstdlib>
#include <iostream>

#ifdef CUDA_MANAGED
#include <cuda_runtime.h>
#endif

#include <cedar/util/basic_allocator.h>

namespace cedar
{

basic_allocator::basic_allocator(global_params & params) :
	allocator(std::malloc), deallocator(std::free)
{
	#ifdef CUDA_MANAGED
	if (params.memory_type == memtype::managed) {
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

}
