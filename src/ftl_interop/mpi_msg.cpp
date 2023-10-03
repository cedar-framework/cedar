#include <ftl/Base.hpp>
#include <ftl/Kernel.hpp>
#include <ftl/Buffer.hpp>
#include <ftl/Cedar.hpp>
#include <ftl/mpi_msg.hpp>

#define MAX_PROCS 27
#define MAX_PATTERNS 6
#define MPI_STATUS_SIZE (sizeof(MPI_Status))

static int msg_version;

template <typename T, int ndim>
struct FortranArrayView {
    T* base;
    int strides[ndim];

    template <typename U>
    void unroll_stride(int i, int running_stride, U&& u) {
        strides[i] = running_stride;
    }

    template <typename U, typename ...Us>
    void unroll_stride(int i, int running_stride, U&& u, Us&&... us) {
        strides[i] = running_stride;
        unroll_stride(i + 1, running_stride * u, us...);
    }

    template <typename ...Us>
    FortranArrayView(T* base, Us&&... us): base(base) {
        unroll_stride(0, 1, us...);
    }

    template <typename U>
    int unroll_idx(int i, U&& u) {
        return u * strides[i];
    }

    template <typename U, typename ...Us>
    int unroll_idx(int i, U&& u, Us&&... us) {
        return u * strides[i] + unroll_idx(i + 1, us...);
    }

    template <typename ...Us>
    T& operator()(Us&&... us) {
        const int idx = unroll_idx(0, us...);
        return base[idx];
    }

    /* Playing fast and loose with the automatic conversions :-) */
    operator MPI_Status&() {
        return *((MPI_Status*) base);
    }
};

extern "C" {
    extern struct {
        int32_t MSG_sendid[MAX_PROCS * MAX_PATTERNS];
        int32_t MSG_recvid[MAX_PROCS * MAX_PATTERNS];
        int32_t SendStatus[MPI_STATUS_SIZE];
        int32_t RecvStatus[MPI_STATUS_SIZE];
        int32_t MSGSegment[MAX_PATTERNS];
        int32_t MSG_COMM;
        int32_t MSG_COMM_PARENT;
        int32_t MSG_COMM_PARENT_FLAG;
        int32_t MSG_TRANSFER_TYPE[MAX_PATTERNS];
        int32_t MSG_BLOCKING;
    } msg_sendrec_;
}

static FortranArrayView<int32_t, 2> msg_recvid(msg_sendrec_.MSG_recvid, MAX_PROCS, MAX_PATTERNS);
static FortranArrayView<int32_t, 2> msg_sendid(msg_sendrec_.MSG_sendid, MAX_PROCS, MAX_PATTERNS);
static FortranArrayView<int32_t, 1> msgsegment(msg_sendrec_.MSGSegment, MAX_PATTERNS);
static FortranArrayView<int32_t, 1> recvstatus(msg_sendrec_.RecvStatus, MPI_STATUS_SIZE);
static FortranArrayView<int32_t, 1> sendstatus(msg_sendrec_.SendStatus, MPI_STATUS_SIZE);
static FortranArrayView<int32_t, 1> msg_transfer_type(msg_sendrec_.MSG_TRANSFER_TYPE, MAX_PATTERNS);
static int32_t& msg_blocking = msg_sendrec_.MSG_BLOCKING;
static int32_t& msg_comm = msg_sendrec_.MSG_COMM;
static int32_t& msg_comm_parent = msg_sendrec_.MSG_COMM_PARENT;
static int32_t& msg_comm_parent_flag = msg_sendrec_.MSG_COMM_PARENT_FLAG;

void MSG_enable(int myproc, int numproc) {
    int ierror;
    int flag;
    msg_version = 2;
    msg_blocking = 1;
    MPI_Initialized(&flag);
    if (!flag) {
        ierror = MPI_Init(nullptr, nullptr);
    }
    if (msg_comm_parent_flag != -1) {
        msg_comm_parent = mpi_comm_world;
        MPI_Comm_dup(mpi_comm_world, msg_comm, ierror);
    } else {
        MPI_Comm_dup(msg_comm_parent, msg_comm, ierror);
    }
    msg_comm_parent_flag = 0;
    MPI_Comm_rank(msg_comm, myproc, ierror);
    myproc ++;
    MPI_Comm_size(msg_comm, numproc, ierror);
}

void MSG_set_comm_parent(int32_t comm) {
    msg_comm_parent_flag = -1;
    msg_comm_parent = comm;
}

int MSG_myproc() {
    int ierror, myproc;
    MPI_Comm_rank(msg_comm, myproc, ierror);
    return myproc + 1;
}

int MSG_nproc() {
    int ierror, numproc;
    MPI_Comm_size(msg_comm, numproc, ierror);
    return numproc;
}

void MSG_comm_type(bool t) {
    msg_blocking = t;
}

void MSG_disable(int& ierror) {

    MPI_Barrier(msg_comm, ierror);
    MPI_Comm_free(msg_comm, ierror);
}

void MSG_tbdx_send(
	ftl::Buffer<real_t> x,
	ftl::Buffer<real_t> y,
	int nproc,
	ftl::Buffer<len_t> proc,
	ftl::Buffer<len_t> ipr,
	ftl::Buffer<len_t> index,
	int ptrn,
	int& ierr) {

    static bool first_call = true;

    int iproc, insegmentstart, insegmentsize, myproc, outsegmentsize, outsegmentstart;
    ierr = 0;
    if (msg_blocking==1) {
        if (nproc==0) {
            return;
        }
        myproc = MSG_myproc();
        if (proc((1) - 1)<0) {
            proc((1) - 1) = (-proc((1) - 1));
        }
        for (iproc = 1; iproc <= nproc; ++iproc) {
            outsegmentsize = (ipr((2*iproc) - 1)-ipr(((2*iproc)-1) - 1));
            outsegmentstart = ipr(((2*iproc)-1) - 1);
            insegmentsize = (ipr(((2*iproc)+1) - 1)-ipr((2*iproc) - 1));
            insegmentstart = ipr((2*iproc) - 1);
            if (myproc<proc((iproc) - 1)) {
                if (index((outsegmentstart) - 1)>=0) {
                    MSG_tbdx_gather(x, y, iproc, ipr, index);
                    MPI_Send(static_cast<void*>(y(0)), outsegmentsize, MPI_DOUBLE, (proc((iproc) - 1)-1), ptrn, msg_comm, ierr);
                    if (ierr!=MPI_SUCCESS) {
                        return;
                    }
                } else {
                    MPI_Send(static_cast<void*>(x((-index((outsegmentstart) - 1)) - 1)), outsegmentsize, MPI_DOUBLE, (proc((iproc) - 1)-1), ptrn, msg_comm, ierr);
                }
                if (index((insegmentstart) - 1)>=0) {
                    MPI_Recv(static_cast<void*>(y(0)), insegmentsize, MPI_DOUBLE, (proc((iproc) - 1)-1), ptrn, msg_comm, recvstatus, ierr);
                    if (ierr!=MPI_SUCCESS) {
                        return;
                    }
                    MSG_tbdx_scatter(x, y, iproc, ipr, index);
                } else {
                    MPI_Recv(static_cast<void*>(x((-index((insegmentstart) - 1)) - 1)), insegmentsize, MPI_DOUBLE, (proc((iproc) - 1)-1), ptrn, msg_comm, recvstatus, ierr);
                    if (ierr!=MPI_SUCCESS) {
                        return;
                    }
                }
            } else {
                if (index((insegmentstart) - 1)>=0) {
                    MPI_Recv(static_cast<void*>(y(0)), insegmentsize, MPI_DOUBLE, (proc((iproc) - 1)-1), ptrn, msg_comm, recvstatus, ierr);
                    if (ierr!=MPI_SUCCESS) {
                        return;
                    }
                    MSG_tbdx_scatter(x, y, iproc, ipr, index);
                } else {
                    MPI_Recv(static_cast<void*>(x((-index((insegmentstart) - 1)) - 1)), insegmentsize, MPI_DOUBLE, (proc((iproc) - 1)-1), ptrn, msg_comm, recvstatus, ierr);
                    if (ierr!=MPI_SUCCESS) {
                        return;
                    }
                }
                if (index((outsegmentstart) - 1)>=0) {
                    MSG_tbdx_gather(x, y, iproc, ipr, index);
                    MPI_Send(static_cast<void*>(y(0)), outsegmentsize, MPI_DOUBLE, (proc((iproc) - 1)-1), ptrn, msg_comm, ierr);
                    if (ierr!=MPI_SUCCESS) {
                        return;
                    }
                } else {
                    MPI_Send(static_cast<void*>(x((-index((outsegmentstart) - 1)) - 1)), outsegmentsize, MPI_DOUBLE, (proc((iproc) - 1)-1), ptrn, msg_comm, ierr);
                }
            }
        }
    } else { /* Use non-blocking communication */
        if (nproc==0) {
            return;
        }
        if ((ptrn>MAX_PATTERNS)||(ptrn<=0)) {
            ierr = -2;
            return;
        }
        if ((nproc>MAX_PROCS)||(nproc<0)) {
            ierr = -3;
            return;
        }
        if (first_call) {
            for (iproc = 1; iproc <= MAX_PATTERNS; ++iproc) {
                msg_sendid(0, iproc - 1) = 0;
                msg_recvid(0, iproc - 1) = 0;
            }
            first_call = false;
        }
        if (msg_sendid(0, ptrn - 1)==0 && msg_recvid(0, ptrn - 1)==0) {
            /*
             Open the communications channels, set the type of data transfer
             for this pattern
             */
            if (proc(0)<0) {
                msg_transfer_type(ptrn - 1) = 1;
                proc(0) = -proc(0);
            } else {
                msg_transfer_type(ptrn - 1) = 0;
            }

            /*
              Find the maximal size of outgoing segment in order to find a
              "safe" place within the array y to put the buffer for the incoming
              data
            */
            msgsegment(ptrn - 1) = 0;
            for (iproc = 1; iproc <= nproc; ++iproc) {
                outsegmentsize = (ipr((2*iproc) - 1)-ipr(((2*iproc)-1) - 1));
                if (msgsegment(ptrn - 1) < outsegmentsize) {
                    msgsegment(ptrn - 1) = outsegmentsize;
                }
            }
            msgsegment(ptrn - 1) = msgsegment(ptrn);

            /* Open up channels */
            for (iproc = 1; iproc <= nproc; ++iproc) {
                outsegmentsize = (ipr((2*iproc) - 1)-ipr(((2*iproc)-1) - 1));
                outsegmentstart = ipr(((2*iproc)-1) - 1);
                if (outsegmentsize>0) {
                    if (index((outsegmentstart) - 1)>=0) {
                        /* Noncontiguous memory segment: give the buffer's address */
                        MPI_Send_init(static_cast<void*>(y(0)), outsegmentsize, MPI_DOUBLE, (proc((iproc) - 1)-1), ptrn, msg_comm, msg_sendid(iproc - 1, ptrn - 1), ierr);
                        if (ierr!=MPI_SUCCESS) {
                            return;
                        }
                    } else {
                        /* Contiguous memory segment: give the data address */
                        MPI_Send_init(static_cast<void*>(x((-index((outsegmentstart) - 1)) - 1)), outsegmentsize, MPI_DOUBLE, (proc((iproc) - 1)-1), ptrn, msg_comm, msg_sendid(iproc - 1, ptrn - 1), ierr);
                        if (ierr!=MPI_SUCCESS) {
                            return;
                        }
                    }
                }

                insegmentsize = (ipr(((2*iproc)+1) - 1)-ipr((2*iproc) - 1));
                insegmentstart = ipr((2*iproc) - 1);
                if (insegmentsize>0) {
                    if (index((insegmentstart) - 1)>=0) {
                        MPI_Recv_init(static_cast<void*>(y((msgsegment(ptrn - 1)) - 1)), insegmentsize, MPI_DOUBLE, (proc((iproc) - 1)-1), ptrn, msg_comm, msg_recvid(iproc - 1, ptrn - 1), ierr);
                        if (ierr!=MPI_SUCCESS) {
                            return;
                        }
                    } else {
                        MPI_Recv_init(static_cast<void*>(x((-index((insegmentstart) - 1)) - 1)), insegmentsize, MPI_DOUBLE, (proc((iproc) - 1)-1), ptrn, msg_comm, msg_recvid(iproc - 1, ptrn - 1), ierr);
                        if (ierr!=MPI_SUCCESS) {
                            return;
                        }
                    }
                }
            }
        }

        if (msg_transfer_type(ptrn - 1)==1) {
            /* Exchange data through the channels using all to all
               Send all messages out one by one */

            for (iproc = 1; iproc <= nproc; ++iproc) {
                outsegmentsize = (ipr((2*iproc) - 1)-ipr(((2*iproc)-1) - 1));
                if (outsegmentsize>0) {
                    /* Gather outgoing data in the outgoing buffer */
                    MSG_tbdx_gather(x, y, iproc, ipr, index);

                    /* Start sending the outgoing data to iproc */
                    MPI_Start(msg_sendid(iproc - 1, ptrn - 1), ierr);
                    if (ierr!=MPI_SUCCESS) {
                        return;
                    }
                    MPI_Wait(msg_sendid(iproc - 1, ptrn - 1), sendstatus, ierr);
                    if (ierr!=MPI_SUCCESS) {
                        return;
                    }
                }
            }
        } else {
            /* Exchange data through the channels using shifts */

            if (nproc>1) {
                for (iproc = 1; iproc <= (nproc-1); ++iproc) {
                    outsegmentsize = (ipr((2*iproc) - 1)-ipr(((2*iproc)-1) - 1));
                    insegmentsize = (ipr(((2*iproc)+1) - 1)-ipr((2*iproc) - 1));
                    /* Start receiving the incoming data from iproc */

                    if (insegmentsize>0) {
                        MPI_Start(msg_recvid(iproc - 1, ptrn - 1), ierr);
                        if (ierr!=MPI_SUCCESS) {
                            return;
                        }
                    }
                    if (outsegmentsize>0) {
                        MSG_tbdx_gather(x, y, iproc, ipr, index);
                        MPI_Start(msg_sendid(iproc - 1, ptrn - 1), ierr);
                        if (ierr!=MPI_SUCCESS) {
                            return;
                        }
                    }
                    if (insegmentsize>0) {
                        MPI_Wait(msg_recvid(iproc - 1, ptrn - 1), recvstatus, ierr);
                        if (ierr!=MPI_SUCCESS) {
                            return;
                        }
                        MSG_tbdx_scatter(x, y(msgsegment(ptrn - 1) - 1), iproc, ipr, index);
                    }
                    if (outsegmentsize>0) {
                        MPI_Wait(msg_sendid(iproc - 1, ptrn - 1), sendstatus, ierr);
                        if (ierr!=MPI_SUCCESS) {
                            return;
                        }
                    }
                }
            }
            iproc = nproc;
            outsegmentsize = (ipr((2*iproc) - 1)-ipr(((2*iproc)-1) - 1));
            insegmentsize = (ipr(((2*iproc)+1) - 1)-ipr((2*iproc) - 1));

            /* Start receiving the incoming data from nproc */
            if (insegmentsize>0) {
                MPI_Start(msg_recvid(iproc - 1, ptrn - 1), ierr);
                if (ierr!=MPI_SUCCESS) {
                    return;
                }
            }
            if (outsegmentsize>0) {
                MSG_tbdx_gather(x, y, iproc, ipr, index);
                MPI_Start(msg_sendid(iproc - 1, ptrn - 1), ierr);
                if (ierr!=MPI_SUCCESS) {
                    return;
                }
            }
        }
    }
}

void MSG_tbdx_close(
	ftl::Buffer<real_t> x,
	ftl::Buffer<real_t> y,
	int nproc,
	ftl::Buffer<len_t> proc,
	ftl::Buffer<len_t> ipr,
	ftl::Buffer<len_t> index,
	int ptrn,
        int& ierr) {

    int ii;
    ierr = 0;
    if (msg_blocking==0) {
        if (nproc==0) {
            return;
        }
        for (ii = 1; ii <= nproc; ++ii) {
            if (msg_sendid(ii - 1, ptrn - 1) != 0) {
                MPI_Request_free(msg_sendid(ii - 1, ptrn - 1), ierr);
                if (ierr!=MPI_SUCCESS) {
                    return;
                }
            }
            if (msg_recvid(ii - 1, ptrn - 1) != 0) {
                MPI_Request_free(msg_recvid(ii - 1, ptrn - 1), ierr);
                if (ierr!=MPI_SUCCESS) {
                    return;
                }
            }
        }
        msg_sendid(0, ptrn - 1) = 0;
        msg_recvid(0, ptrn - 1) = 0;
        ierr = 0;
    }
}

void MSG_tbdx_receive(
	ftl::Buffer<real_t> x,
	ftl::Buffer<real_t> y,
	int nproc,
	ftl::Buffer<len_t> proc,
	ftl::Buffer<len_t> ipr,
	ftl::Buffer<len_t> index,
	int ptrn,
	int& ierr) {

    int outsegmentsize, insegmentsize, iproc;
    ierr = 0;
    if (msg_blocking==0) {
        if (nproc==0) {
            return;
        }
        if (msg_transfer_type(ptrn - 1) == 1) {
            for (iproc = 1; iproc <= nproc; ++iproc) {
                insegmentsize = (ipr(((2*iproc)+1) - 1)-ipr((2*iproc) - 1));
                if (insegmentsize>0) {
                    MPI_Start(msg_recvid(iproc - 1, ptrn - 1), ierr);
                    MPI_Wait(msg_recvid(iproc - 1, ptrn - 1), recvstatus, ierr);
                    MSG_tbdx_scatter(x, y(msgsegment(ptrn - 1) - 1), iproc, ipr, index);
                }
            }
        } else {
            outsegmentsize = (ipr((2*nproc) - 1)-ipr(((2*nproc)-1) - 1));
            insegmentsize = (ipr(((2*nproc)+1) - 1)-ipr((2*nproc) - 1));
            if (insegmentsize>0) {
                MPI_Wait(msg_recvid(nproc - 1, ptrn - 1), recvstatus, ierr);
                MSG_tbdx_scatter(x, y(msgsegment(ptrn - 1) - 1), nproc, ipr, index);
            }
            if (outsegmentsize>0) {
                MPI_Wait(msg_sendid(nproc - 1, ptrn - 1), sendstatus, ierr);
                if (ierr!=MPI_SUCCESS) {
                    return;
                }
            }
        }
    }
}

void MSG_tbdx_gather(
	ftl::Buffer<real_t> x,
	ftl::Buffer<real_t> y,
	int iproc,
	ftl::Buffer<len_t> ipr,
	ftl::Buffer<len_t> index) {

    int outsegmentstart, outsegmentend, outsegmentsize, j, k;
    outsegmentstart = ipr(((2*iproc)-1) - 1);
    if (index((outsegmentstart) - 1)<0) {
        return;
    }
    outsegmentsize = (ipr((2*iproc) - 1)-ipr(((2*iproc)-1) - 1));
    outsegmentend = ((outsegmentstart+outsegmentsize)-1);
    k = 1;
    for (j = outsegmentstart; j <= outsegmentend; ++j) {
        y((k) - 1) = x((index((j) - 1)) - 1);
        k = (k+1);
    }
}

void MSG_tbdx_scatter(
	ftl::Buffer<real_t> x,
	ftl::Buffer<real_t> y,
	int iproc,
	ftl::Buffer<len_t> ipr,
	ftl::Buffer<len_t> index) {

    int insegmentsize, insegmentstart, insegmentend, j, k;
    insegmentstart = ipr((2*iproc) - 1);
    if (index((insegmentstart) - 1)<0) {
        return;
    }
    insegmentsize = (ipr(((2*iproc)+1) - 1)-ipr((2*iproc) - 1));
    insegmentend = ((insegmentstart+insegmentsize)-1);
    k = 1;
    for (j = insegmentstart; j <= insegmentend; ++j) {
        x((index((j) - 1)) - 1) = y((k) - 1);
        k = (k+1);
    }
}
