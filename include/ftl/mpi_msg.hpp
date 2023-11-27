#ifndef FTL_CONVERTED_MPI_MSG_INC_
#define FTL_CONVERTED_MPI_MSG_INC_

#include <ftl/Base.hpp>
#include <ftl/Kernel.hpp>
#include <ftl/Buffer.hpp>
#include <ftl/Cedar.hpp>

#include <mpi.h>

void MSG_enable(
	int myproc,
	int numproc);

void MSG_set_comm_parent(int32_t comm);

int MSG_myproc();

int MSG_nproc();

void MSG_comm_type(bool t);

void MSG_disable(int& ierror);

void MSG_tbdx_send(
	ftl::Buffer<real_t> x,
	ftl::Buffer<real_t> y,
	int nproc,
	ftl::Buffer<len_t> proc,
	ftl::Buffer<len_t> ipr,
	ftl::Buffer<len_t> index,
	int ptrn,
	int& ierr);

void MSG_tbdx_close(
	ftl::Buffer<real_t> x,
	ftl::Buffer<real_t> y,
	int nproc,
	ftl::Buffer<len_t> proc,
	ftl::Buffer<len_t> ipr,
	ftl::Buffer<len_t> index,
	int ptrn,
	int& ierr);

void MSG_tbdx_receive(
	ftl::Buffer<real_t> x,
	ftl::Buffer<real_t> y,
	int nproc,
	ftl::Buffer<len_t> proc,
	ftl::Buffer<len_t> ipr,
	ftl::Buffer<len_t> index,
	int ptrn,
	int& ierr);

void MSG_tbdx_gather(
	ftl::Buffer<real_t> x,
	ftl::Buffer<real_t> y,
	int iproc,
	ftl::Buffer<len_t> ipr,
	ftl::Buffer<len_t> index);

void MSG_tbdx_scatter(
	ftl::Buffer<real_t> x,
	ftl::Buffer<real_t> y,
	int iproc,
	ftl::Buffer<len_t> ipr,
	ftl::Buffer<len_t> index);

#endif /* FTL_CONVERTED_MPI_MSG_INC_*/
