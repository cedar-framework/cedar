#ifndef FTL_CEDAR_HPP_INC_
#define FTL_CEDAR_HPP_INC_

/* FTL <-> Cedar interop */

#include <cedar/types.h>
#include <cedar/services/halo_exchange.h>
#include <cedar/services/mempool.h>
#include <cstdint>
#include <iostream>

using real_t = cedar::real_t;
using len_t = cedar::len_t;

#define real real_t

#include <ftl/Buffer.hpp>
#include <ftl/mpi.hpp>

#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_PtrMSGSO.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_PtrMSG.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_MSGGrid.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_MSGGridSO.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_PtrLS.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_LSGrid.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_cg_boxmg.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_COPY_cg_WS_SO.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_COPY_cg_WS_RHS.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_COPY_cg_rV_G_L.f90.hpp>
#include <src/2d/ftn/mpi/comm_coarsen.f90.hpp>

template <typename int_type>
void dpttrf(int N,
            ftl::BufferView<real_t> D,
            ftl::BufferView<real_t> E,
            int_type& info) {
    throw std::runtime_error("DPTTRF not implemented");
}
#define DPTTRF dpttrf

template <typename int_type>
void dpttrs(int N,
            int NRHS,
            ftl::BufferView<real_t> D,
            ftl::BufferView<real_t> E,
            ftl::BufferView<real_t> B,
            int ldb,
            int_type& info) {
    throw std::runtime_error("DPTTRF not implemented");
}
#define DPTTRS dpttrs

template <typename int_type>
void dpbtrs(const char* uplo,
            int N,
            int KD,
            int NRHS,
            ftl::BufferView<real_t> AB,
            int LDAB,
            ftl::BufferView<real_t> B,
            int LDB,
            int_type& info) {

    throw std::runtime_error("DPBTRS not implemented");
}
#define DPBTRS dpbtrs

template <typename int_type>
void dpbtrf(const char* uplo,
            int N,
            int KD,
            ftl::BufferView<real_t> AB,
            int LDAB,
            int_type& info) {

    throw std::runtime_error("DPBTRF not implemented");
}
#define DPBTRF dpbtrf

void halo_exchange(int k, ftl::BufferView<real_t> q, void* halof_void) {
    auto* halof =
        static_cast<cedar::services::halo_exchange_base*>(halof_void);
    q->dev_to_host();
    halof->exchange_func(k, q->data());
}

void halo_stencil_exchange(int k, ftl::BufferView<real_t> soc, void *halof_void) {
    auto* halof =
        static_cast<cedar::services::halo_exchange_base*>(halof_void);
    soc->dev_to_host();
    halof->exchange_sten(k, soc->data());
}

void MSG_tp_setup(void* la_size, void* eid_s, void* gc_ld, void* gc_eid,
                  int numproc, int myproc, int nproc, void* proc, void* ipr, void* index,
                  int sfa, int pfa, int& ier) {
    throw std::runtime_error("MSG_tp_setup not implemented");
}

void cedar_mempool_pos(void* mp_void, int nbytes, int& pos) {
    auto* mp =
        static_cast<cedar::services::mempool*>(mp_void);
    std::size_t nbytes_in = nbytes;
    pos = mp->pos(nbytes_in);
}

void print_error(const char* str) {
    std::cout << str << std::endl;
    std::exit(1);
}

#endif
