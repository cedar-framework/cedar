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


void basic_allocator::prefetch_impl(void *addr, std::size_t nbytes)
{
	#ifdef CUDA_MANAGED
	int defdev = omp_get_default_device();
	cudaMemPrefetchAsync(addr, nbytes, defdev, cudaStreamLegacy);
	#endif
}


void basic_allocator::sync()
{
	#ifdef CUDA_MANAGED
	cudaDeviceSynchronize();
	#endif
}

}
