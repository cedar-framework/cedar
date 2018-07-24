subroutine test_gather(tvals, plane_len, recv, comm, mp) BIND(C, name="test_gather")
  use ModInterface
  use message_passing
  implicit none
  include "mpif.h"

  real(real_t) :: tvals(plane_len), recv(plane_len)
  integer(c_int), value :: plane_len
  integer, value :: comm
  type(c_ptr) :: mp
  integer :: ierr

  call cedar_gather(mp, tvals, plane_len, MPI_DOUBLE_PRECISION,&
       recv, plane_len, MPI_DOUBLE_PRECISION, 0, comm, ierr)

end subroutine test_gather
