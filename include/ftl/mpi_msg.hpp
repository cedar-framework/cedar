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

void MSG_set_comm_parent(MPI_Comm comm);

void MSG_set_comm_parent(std::intptr_t comm);

int MSG_myproc();

int MSG_nproc();

void MSG_comm_type(bool t);

void MSG_disable(int& ierror);

void MSG_tbdx_send(
	ftl::BufferView<real_t> x,
	ftl::BufferView<real_t> y,
	int nproc,
	ftl::BufferView<len_t> proc,
	ftl::BufferView<len_t> ipr,
	ftl::BufferView<len_t> index,
	int ptrn,
	int ierr);

void MSG_tbdx_close(
	ftl::BufferView<real_t> x,
	ftl::BufferView<real_t> y,
	int nproc,
	ftl::BufferView<len_t> proc,
	ftl::BufferView<len_t> ipr,
	ftl::BufferView<len_t> index,
	int ptrn,
	int ierr);

void MSG_tbdx_receive(
	ftl::BufferView<real_t> x,
	ftl::BufferView<real_t> y,
	int nproc,
	ftl::BufferView<len_t> proc,
	ftl::BufferView<len_t> ipr,
	ftl::BufferView<len_t> index,
	int ptrn,
	int ierr);

void MSG_tbdx_gather(
	ftl::BufferView<real_t> x,
	ftl::BufferView<real_t> y,
	int iproc,
	ftl::BufferView<len_t> ipr,
	ftl::BufferView<len_t> index);

void MSG_tbdx_scatter(
	ftl::BufferView<real_t> x,
	ftl::BufferView<real_t> y,
	int iproc,
	ftl::BufferView<len_t> ipr,
	ftl::BufferView<len_t> index);

#endif /* FTL_CONVERTED_MPI_MSG_INC_*/
