#include <cedar/services/message_passing.h>

#include <ftl/Cedar.hpp>
#include <unordered_map>

void MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source,
              int tag, int32_t comm, MPI_Status& status, int32_t& ierr) {

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Recv(buf, count, datatype, source, tag, comm_ptr, &status);
}

void MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source,
              int tag, int32_t comm, void* status, int32_t& ierr) {

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Recv(buf, count, datatype, source, tag, comm_ptr, static_cast<MPI_Status*>(status));
}

void MPI_Recv_init(void* buf, int count, MPI_Datatype datatype, int source,
                   int tag, int32_t comm, MPI_Request& request, int32_t& ierr) {

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Recv_init(buf, count, datatype, source, tag, comm_ptr, &request);
}

void MPI_Recv_init(void* buf, int count, MPI_Datatype datatype, int source,
                   int tag, int32_t comm, int32_t& request, int32_t& ierr) {

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    MPI_Request request_ptr = MPI_Request_f2c(request);
    ierr = MPI_Recv_init(buf, count, datatype, source, tag, comm_ptr, &request_ptr);
    request = MPI_Request_c2f(request_ptr);
}

void MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest,
              int tag, int32_t comm, int32_t& ierr) {

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Send(buf, count, datatype, dest, tag, comm_ptr);
}

void MPI_Send_init(void* buf, int count, MPI_Datatype datatype, int dest,
                   int tag, int32_t comm, MPI_Request& request, int32_t& ierr) {

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Send_init(buf, count, datatype, dest, tag, comm_ptr, &request);
}

void MPI_Send_init(void* buf, int count, MPI_Datatype datatype, int dest,
                   int tag, int32_t comm, int32_t& request, int32_t& ierr) {

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    MPI_Request request_ptr = MPI_Request_f2c(request);
    ierr = MPI_Send_init(buf, count, datatype, dest, tag, comm_ptr, &request_ptr);
    request = MPI_Request_c2f(request_ptr);
}

void MPI_Wait(MPI_Request& request, MPI_Status& status, int32_t& ierr) {
    ierr = MPI_Wait(&request, &status);
}

void MPI_Wait(MPI_Request& request, void* status, int32_t& ierr) {
    ierr = MPI_Wait(&request, static_cast<MPI_Status*>(status));
}

void MPI_Wait(int32_t& request, MPI_Status& status, int32_t& ierr) {
    MPI_Request request_ptr = MPI_Request_f2c(request);
    ierr = MPI_Wait(&request_ptr, &status);
    request = MPI_Request_c2f(request_ptr);
}

void MPI_Wait(int32_t& request, void* status, int32_t& ierr) {
    MPI_Request request_ptr = MPI_Request_f2c(request);
    ierr = MPI_Wait(&request_ptr, static_cast<MPI_Status*>(status));
    request = MPI_Request_c2f(request_ptr);
}

void MPI_Sendrecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                  void* recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                  int32_t comm, MPI_Status& status, int32_t& ierr) {

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag,
                        recvbuf, recvcount, recvtype, source, recvtag,
                        comm_ptr, &status);

}

void MPI_Sendrecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
                  void* recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
                  int32_t comm, void* status, int32_t& ierr) {

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag,
                        recvbuf, recvcount, recvtype, source, recvtag,
                        comm_ptr, static_cast<MPI_Status*>(status));

}

void MPI_Allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                   void* recvbuf, int recvcount, MPI_Datatype recvtype,
                   int32_t comm, int32_t& ierr) {

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Allgather(sendbuf, sendcount, sendtype,
                         recvbuf, recvcount, recvtype,
                         comm_ptr);
}

void MPI_Allgather(ftl::Buffer<real_t> sendbuf, int sendcount, MPI_Datatype sendtype,
                   ftl::Buffer<real_t> recvbuf, int recvcount, MPI_Datatype recvtype,
                   int32_t comm, int32_t& ierr) {

    sendbuf.dev_to_host();
    recvbuf.mark_device_dirty(false);
    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);

    ierr = MPI_Allgather(sendbuf.data(), sendcount, sendtype,
                         recvbuf.data(), recvcount, recvtype,
                         comm_ptr);

    recvbuf.mark_host_dirty(true);
}

void MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                 void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,
                 int32_t comm, int32_t& ierr) {

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Scatter(sendbuf, sendcount, sendtype,
                       recvbuf, recvcount, recvtype,
                       root, comm_ptr);

}

void MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,
                int32_t comm, int32_t& ierr) {

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Gather(sendbuf, sendcount, sendtype,
                      recvbuf, recvcount, recvtype,
                      root, comm_ptr);

}

void MPI_Gather(ftl::Buffer<real_t> sendbuf, int sendcount, MPI_Datatype sendtype,
                ftl::Buffer<real_t> recvbuf, int recvcount, MPI_Datatype recvtype, int root,
                int32_t comm, int32_t& ierr) {

    sendbuf.dev_to_host();
    recvbuf.mark_device_dirty(false);

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Gather(sendbuf.data(), sendcount, sendtype,
                      recvbuf.data(), recvcount, recvtype,
                      root, comm_ptr);

    recvbuf.mark_host_dirty(true);
}

void MPI_Comm_rank(int32_t comm, int32_t& rank, int32_t& ierr) {
    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Comm_rank(comm_ptr, &rank);
}

void MPI_Comm_dup(int32_t comm, int32_t& new_comm, int32_t& ierr) {
    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    MPI_Comm new_comm_ptr;
    ierr = MPI_Comm_dup(comm_ptr, &new_comm_ptr);
    new_comm = MPI_Comm_c2f(new_comm_ptr);
}

void MPI_Comm_free(int32_t comm, int32_t& ierr) {
    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Comm_free(&comm_ptr);
}

void MPI_Comm_size(int32_t comm, int32_t& size, int32_t& ierr) {
    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Comm_size(comm_ptr, &size);
}

void MPI_Request_free(MPI_Request& request, int32_t& ierr) {
    ierr = MPI_Request_free(&request);
}

void MPI_Request_free(int32_t& request, int32_t& ierr) {
    MPI_Request request_ptr = MPI_Request_f2c(request);
    ierr = MPI_Request_free(&request_ptr);
    request = MPI_Request_c2f(request_ptr);
}

void MPI_Start(MPI_Request& request, int32_t& ierr) {
    ierr = MPI_Start(&request);
}

void MPI_Start(int32_t& request, int32_t& ierr) {
    MPI_Request request_ptr = MPI_Request_f2c(request);
    ierr = MPI_Start(&request_ptr);
    request = MPI_Request_c2f(request_ptr);
}

void MPI_Barrier(int32_t comm, int32_t& ierr) {
    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Barrier(comm_ptr);
}

void MPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
               int root, int32_t comm, int32_t& ierr) {

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Bcast(buffer, count, datatype, root, comm_ptr);
}

void MPI_Bcast(ftl::Buffer<real_t> buffer, int count, MPI_Datatype datatype,
               int root, int32_t comm, int32_t& ierr) {

    buffer.dev_to_host();

    MPI_Comm comm_ptr = MPI_Comm_f2c(comm);
    ierr = MPI_Bcast(buffer.data(), count, datatype, root, comm_ptr);
}

void MPI_Finalize(int32_t comm) {
    MPI_Finalize();
}

void cedar_comm_split(void* mp_void, int fcomm, int color, int key, int& fnewcomm, int& ierr) {
    auto mp = static_cast<cedar::services::message_passing*>(mp_void);

    MPI_Comm comm = MPI_Comm_f2c(fcomm);
    MPI_Comm newcomm;
    ierr = mp->comm_split(comm, color, key, &newcomm);
    fnewcomm = MPI_Comm_c2f(newcomm);
}
