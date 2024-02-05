#include <ftl/Cedar.hpp>
#include <ftl/Blas.hpp>
#include <iostream>
#include <cedar/services/mempool.h>

void halo_exchange(int k, ftl::Buffer<real_t> q, void* halof_void) {
    auto halof = static_cast<cedar::services::halo_exchange_base*>(halof_void);
    halof->exchange_func(k, q);
}

void halo_stencil_exchange(int k, ftl::Buffer<real_t> soc, void *halof_void) {
    auto halof = static_cast<cedar::services::halo_exchange_base*>(halof_void);
    halof->exchange_sten(k, soc);
}

void cedar_mempool_pos(void* mp_void, int nbytes, int& pos) {
    auto mp = static_cast<cedar::services::mempool*>(mp_void);
    pos = mp->pos(nbytes);
}

void dpotrf_gpu(const char* uplo, int n, ftl::Buffer<real_t>& abd, int nabd, int& info) {
    abd.host_to_dev();

    ftl::blas::TriangularFillMode fill_mode =
        (uplo[0] == 'U' ?
         ftl::blas::TriangularFillMode::Upper :
         ftl::blas::TriangularFillMode::Lower);

    ftl::blas::cholesky_factorize(fill_mode, abd);
    info = 0;

    abd.mark_device_dirty();
}

void dpotrs_gpu(const char* uplo, int n, int nrhs, ftl::Buffer<real_t>& A, int lda,
                ftl::FlatBuffer<real_t>& B, int ldb, int& info) {
    A.host_to_dev();
    B.host_to_dev();

    ftl::blas::TriangularFillMode fill_mode =
        (uplo[0] == 'U' ?
         ftl::blas::TriangularFillMode::Upper :
         ftl::blas::TriangularFillMode::Lower);

    ftl::blas::cholesky_solve(fill_mode, A, nrhs, B);
    info = 0;

    B.mark_device_dirty();
}

void print_error(const char* str) {
    std::cerr << str << std::endl;
    std::exit(1);
}
