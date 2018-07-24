module message_passing

  use iso_c_binding, only: c_int, c_ptr
  interface
     subroutine cedar_comm_split(mp, comm, color, key, newcomm, ierr) &
          bind(C, name="cedar_comm_split")
       use iso_c_binding, only: c_int, c_ptr
       integer(c_int), value :: color, key
       integer, value :: comm
       integer :: newcomm, ierr
       type(c_ptr) :: mp
     end subroutine cedar_comm_split
  end interface

  interface
     subroutine cedar_gather_impl(mp, sendbuf, sendcount, sendtype,&
          recvbuf, recvcount, recvtype,&
          root, comm, ierr) BIND(C, name="cedar_gather_impl")
       use iso_c_binding, only: c_int, c_ptr
       integer, value :: comm, sendtype, recvtype
       integer(c_int), value :: sendcount, recvcount, root
       type(c_ptr), value :: sendbuf, recvbuf
       integer :: ierr
       type(c_ptr) :: mp
     end subroutine cedar_gather_impl
  end interface

  interface
     subroutine cedar_scatter_impl(mp, sendbuf, sendcount, sendtype,&
          recvbuf, recvcount, recvtype,&
          root, comm, ierr) BIND(C, name="cedar_scatter_impl")
       use iso_c_binding, only: c_int, c_ptr
       integer, value :: comm, sendtype, recvtype
       integer(c_int), value :: sendcount, recvcount, root
       type(c_ptr), value :: sendbuf, recvbuf
       integer :: ierr
       type(c_ptr) :: mp
     end subroutine cedar_scatter_impl
  end interface

contains
  subroutine cedar_gather(mp, sendbuf, sendcount, sendtype,&
       recvbuf, recvcount, recvtype,&
       root, comm, ierr)
    use ModInterface
    use iso_c_binding, only : c_int, c_ptr, c_loc
    integer, value :: comm, sendtype, recvtype
    integer(c_int), value :: sendcount, recvcount, root
    real(real_t),target :: sendbuf(sendcount), recvbuf(recvcount)
    integer :: ierr
    type(c_ptr) :: mp

    call cedar_gather_impl(mp, c_loc(sendbuf(1)), sendcount, sendtype,&
         c_loc(recvbuf(1)), recvcount, recvtype, root, comm, ierr)

  end subroutine cedar_gather


  subroutine cedar_scatter(mp, sendbuf, sendcount, sendtype,&
       recvbuf, recvcount, recvtype,&
       root, comm, ierr)
    use ModInterface
    use iso_c_binding, only : c_int, c_ptr, c_loc
    integer, value :: comm, sendtype, recvtype
    integer(c_int), value :: sendcount, recvcount, root
    real(real_t),target :: sendbuf(sendcount), recvbuf(recvcount)
    integer :: ierr
    type(c_ptr) :: mp

    call cedar_scatter_impl(mp, c_loc(sendbuf(1)), sendcount, sendtype,&
         c_loc(recvbuf(1)), recvcount, recvtype, root, comm, ierr)

  end subroutine cedar_scatter

end module message_passing
