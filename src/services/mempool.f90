module mempool

  use iso_c_binding, only: c_int, c_ptr
  interface
     subroutine cedar_mempool_pos(mp, nbytes, pos) &
          bind(C, name="cedar_mempool_pos")
       use iso_c_binding, only: c_int, c_ptr
       integer(c_int), value :: nbytes
       type(c_ptr) :: mp
       integer(c_int) :: pos
     end subroutine cedar_mempool_pos
  end interface

end module mempool
