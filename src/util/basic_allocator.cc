#include <cstdlib>
#include <iostream>

#ifdef OFFLOAD
#include <cuda_runtime.h>
#include <omp.h>
#endif

#include <cedar/util/basic_allocator.h>

namespace cedar
{

basic_allocator::basic_allocator(global_params & params) :
	pfsize(0), allocator(std::malloc), deallocator(std::free)
{
	#ifdef OFFLOAD
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
	#ifdef OFFLOAD
	int defdev = omp_get_default_device();
	cudaMemPrefetchAsync(addr, nbytes, defdev, cudaStreamLegacy);
	hints[addr] = memory::location::gpu;
	pfsize += nbytes;
	#endif
}


void basic_allocator::sync()
{
	#ifdef OFFLOAD
	cudaDeviceSynchronize();
	#endif
}


memory::location basic_allocator::hint_impl(const void *addr)
{
	#ifdef OFFLOAD
	auto it = hints.find(addr);
	if (it != hints.end())
		return it->second;
	else
		return memory::location::cpu;
	#else
	return memory::location::cpu;
	#endif
}

}
