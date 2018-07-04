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

end module message_passing
