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
            ftl::Buffer<real_t> D,
            ftl::Buffer<real_t> E,
            int_type& info) {
    throw std::runtime_error("DPTTRF not implemented");
}
#define DPTTRF dpttrf

template <typename int_type>
void dpttrs(int N,
            int NRHS,
            ftl::Buffer<real_t> D,
            ftl::Buffer<real_t> E,
            ftl::Buffer<real_t> B,
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
            ftl::Buffer<real_t> AB,
            int LDAB,
            ftl::Buffer<real_t> B,
            int LDB,
            int_type& info) {

    throw std::runtime_error("DPBTRS not implemented");
}
#define DPBTRS dpbtrs

template <typename int_type>
void dpbtrf(const char* uplo,
            int N,
            int KD,
            ftl::Buffer<real_t> AB,
            int LDAB,
            int_type& info) {

    throw std::runtime_error("DPBTRF not implemented");
}
#define DPBTRF dpbtrf

template <typename ...Us>
void noop(Us&&... args) {}

void halo_exchange(int k, ftl::Buffer<real_t> q, void* halof_void);
void halo_stencil_exchange(int k, ftl::Buffer<real_t> soc, void *halof_void);
void MSG_tp_setup(ftl::Buffer<len_t> la_size, ftl::Buffer<len_t> eid_s, ftl::Buffer<len_t> gc_ld, ftl::Buffer<len_t> gc_eid,
                  int numproc, int myproc, int nproc, ftl::Buffer<len_t> proc, ftl::Buffer<len_t> ipr, ftl::Buffer<len_t> index,
                  int sfa, int pfa, int& ier);
void cedar_mempool_pos(void* mp_void, int nbytes, int& pos);
void print_error(const char* str);

#endif
