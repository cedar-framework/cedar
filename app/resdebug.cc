#include <time.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include <omp.h>
#include <cuda_runtime.h>

#include <cedar/2d/types.h>

using namespace cedar;
using namespace cedar::cdr2;
typedef double real;

#define N 1500000
#define NRUNS 10
#define MANAGED

static grid_func allocvec(len_t nx, len_t ny)
{
	real *addr;
	cudaMallocManaged((void**)&addr, (nx+2)*(ny+2)*sizeof(real), cudaMemAttachGlobal);
	grid_func ret(addr, nx, ny);

	return ret;
}

static void prefetch(grid_func &x)
{
	int defdev = omp_get_default_device();
	cudaMemPrefetchAsync(x.data(), x.size()*sizeof(real), defdev, cudaStreamLegacy);
}


void set_constant(grid_func & x, real val)
{
	int i;
	real *data = x.data();
	#pragma omp target teams distribute parallel for simd
	for (i = 0; i < x.size(); ++i) {
		data[i] = val;
	}

}

int main(int argc, char *argv[])
{
	int defdev = omp_get_default_device();
	real *A;

	config conf("restest.json");
	auto ndofs = conf.getvec<len_t>("grid.n");
	len_t nx = 1200;
	len_t ny = 1200;

	auto x = allocvec(nx, ny);
	struct timespec beg, end;

	int i, j;
	prefetch(x);
	cudaDeviceSynchronize();

	set_constant(x, 1);
	// initialize data

	return 0;
}
