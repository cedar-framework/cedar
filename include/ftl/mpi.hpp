#ifndef FTL_MPI_HPP_INC_
#define FTL_MPI_HPP_INC_

#include <ftl/mpi_msg.hpp>

#define mpi_double_precision MPI_DOUBLE_PRECISION
#define mpi_in_place MPI_IN_PLACE
#define mpi_undefined MPI_UNDEFINED
#define mpi_comm_null (MPI_Comm_c2f(MPI_COMM_NULL))
#define mpi_comm_world (MPI_Comm_c2f(MPI_COMM_WORLD))
#define mpi_comm_self (MPI_Comm_c2f(MPI_COMM_SELF))
#define mpi_status_size (sizeof(MPI_Status) / sizeof(int32_t))

void MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source,
              int tag, int32_t comm, MPI_Status& status, int32_t& ierr);

void MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source,
              int tag, int32_t comm, void* status, int32_t& ierr);

void MPI_Recv_init(void* buf, int count, MPI_Datatype datatype, int source,
                   int tag, int32_t comm, MPI_Request& request, int32_t& ierr);
void MPI_Recv_init(void* buf, int count, MPI_Datatype datatype, int source,
                   int tag, int32_t comm, int32_t& request, int32_t& ierr);

void MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest,
              int tag, int32_t comm, int32_t& ierr);

void MPI_Send_init(void* buf, int count, MPI_Datatype datatype, int dest,
                   int tag, int32_t comm, MPI_Request& request, int32_t& ierr);
void MPI_Send_init(void* buf, int count, MPI_Datatype datatype, int dest,
                   int tag, int32_t comm, int32_t& request, int32_t& ierr);

void MPI_Wait(MPI_Request& request, MPI_Status& status, int32_t& ierr);
void MPI_Wait(MPI_Request& request, void* status, int32_t& ierr);
void MPI_Wait(int32_t& request, MPI_Status& status, int32_t& ierr);
void MPI_Wait(int32_t& request, void* status, int32_t& ierr);

void MPI_Sendrecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                  void* recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                  int32_t comm, MPI_Status& status, int32_t& ierr);

void MPI_Sendrecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                  void* recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                  int32_t comm, void* status, int32_t& ierr);

void MPI_Allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                   void* recvbuf, int recvcount, MPI_Datatype recvtype,
                   int32_t comm, int32_t& ierr);

void MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                 void* recvbuf, int recvcount, MPI_Datatype recvtype,
                 int root, int32_t comm, int32_t& ierr);

void MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                void* recvbuf, int recvcount, MPI_Datatype recvtype,
                int root, int32_t comm, int32_t& ierr);

void MPI_Comm_rank(int32_t comm, int32_t& rank, int32_t& ierr);

void MPI_Comm_dup(int32_t comm, int32_t& new_comm, int32_t& ierr);

void MPI_Comm_free(int32_t comm, int32_t& ierr);

void MPI_Comm_size(int32_t comm, int32_t& size, int32_t& ierr);

void MPI_Request_free(MPI_Request& request, int32_t& ierr);
void MPI_Request_free(int32_t& request, int32_t& ierr);

void MPI_Start(MPI_Request& request, int32_t& ierr);
void MPI_Start(int32_t& request, int32_t& ierr);

void MPI_Barrier(int32_t comm, int32_t& ierr);

void MPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
               int root, int32_t comm, int32_t& ierr);

void MPI_Finalize(int32_t comm);

void cedar_comm_split(void* mp_void, int fcomm, int color, int key, int& fnewcomm, int& ierr);

#endif
