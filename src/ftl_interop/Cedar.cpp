#include <ftl/Cedar.hpp>
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

void print_error(const char* str) {
    std::cerr << str << std::endl;
    std::exit(1);
}
