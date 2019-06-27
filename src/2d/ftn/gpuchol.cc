#include <cedar/types.h>
#include <cedar/memory.h>
#include <cedar/log.h>
#include <cusolverDn.h>
#include <cuda_runtime.h>

extern "C" {

	void init_lapack_gpu()
	{
		using namespace cedar;
		auto *handle = memory::alloc<cusolverDnHandle_t>(1);
		cusolverDnCreate(handle);
		memory::save<cusolverDnHandle_t>("cusolverDnHandle", handle);
	}

	void dpotrs_gpu(char uplo, int N, int NRHS, double *A, int lda, double *B, int ldb, int *info)
	{
		using namespace cedar;
		auto *handle = memory::load<cusolverDnHandle_t>("cusolverDnHandle");
		if (handle == nullptr) {
			log::error << "lapack_gpu not intialized!" << std::endl;
			return;
		}
		cublasFillMode_t fillmode = CUBLAS_FILL_MODE_UPPER;
		if (uplo == 'L')
			fillmode = CUBLAS_FILL_MODE_LOWER;
		auto status = cusolverDnDpotrs(*handle, fillmode, N, NRHS, A, lda, B, ldb, info);
		if (status != CUSOLVER_STATUS_SUCCESS)
			log::error << "Coarse cholesky solve error." << std::endl;
		cudaDeviceSynchronize();
	}

	void dpotrf_gpu(char uplo, int N, double *A, int lda, int *info)
	{
		using namespace cedar;
		auto *handle = memory::load<cusolverDnHandle_t>("cusolverDnHandle");
		if (handle == nullptr) {
			log::error << "lapack_gpu not intialized!" << std::endl;
			return;
		}
		cublasFillMode_t fillmode = CUBLAS_FILL_MODE_UPPER;
		if (uplo == 'L')
			fillmode = CUBLAS_FILL_MODE_LOWER;
		int lwork;
		cusolverDnDpotrf_bufferSize(*handle, fillmode, N, A, lda, &lwork);
		auto * workspace = memory::alloc<double>(lwork);
		cusolverDnDpotrf(*handle, fillmode, N, A, lda, workspace, lwork, info);
		memory::free(workspace);
	}

	void finalize_lapack_gpu()
	{
		using namespace cedar;
		auto *handle = memory::load<cusolverDnHandle_t>("cusolverDnHandle");
		if (handle != nullptr) {
			cusolverDnDestroy(*handle);
			memory::free(handle);
		}
	}
}
