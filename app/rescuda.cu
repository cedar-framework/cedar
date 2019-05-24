#include <cuda.h>

#include <cedar/capi.h>

__global__ void res_kern9(int ilen, int jlen, double *r, const double *so, const double *x, const double *b)
{
		const int idx = blockDim.x * blockIdx.x + threadIdx.x;
		const int i = (idx % (ilen-2)) + 1; // start at index 1 (outside fictitious points)
		const int j = (idx / (ilen-2)) + 1; // start at index 1 (outside fictitious points)
		r[j*ilen + i] = b[j*ilen + i]
				+ so[BMG2_W *(ilen*jlen) + j*ilen + i]         * x[j    *ilen + (i-1)]
				+ so[BMG2_W *(ilen*jlen) + j*ilen + (i+1)]     * x[j    *ilen + (i+1)]
				+ so[BMG2_S *(ilen*jlen) + j*ilen + i]         * x[(j-1)*ilen + i]
				+ so[BMG2_S *(ilen*jlen) + (j+1)*ilen + i]     * x[(j+1)*ilen + i]
				+ so[BMG2_SW*(ilen*jlen) + j*ilen + i]         * x[(j-1)*ilen + (i-1)]
				+ so[BMG2_NW*(ilen*jlen) + j*ilen + (i+1)]     * x[(j-1)*ilen + (i+1)]
				+ so[BMG2_NW*(ilen*jlen) + (j+1)*ilen + i]     * x[(j+1)*ilen + (i-1)]
				+ so[BMG2_SW*(ilen*jlen) + (j+1)*ilen + (i+1)] * x[(j+1)*ilen + (i+1)]
				+ so[BMG2_C *(ilen*jlen) + j*ilen + i]         * x[j*ilen + i];
}


__global__ void init_constant(double *arr, double val)
{
		const int idx = blockDim.x * blockIdx.x + threadIdx.x;
		arr[idx] = val;
}


void residual_cuda9(int ilen, int jlen, double *r, const double *so, const double *x, const double *b)
{
		int interior_size = (jlen-2) * (ilen-2);
		res_kern9<<<interior_size/1024, 1024>>>(ilen, jlen, r, so, x, b);
}


void set_constant_cuda(double *arr, int arrsize, double val)
{
		init_constant<<<arrsize/1024, 1024>>>(arr, val);
}
