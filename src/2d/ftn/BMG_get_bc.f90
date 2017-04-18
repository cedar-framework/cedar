subroutine BMG_get_bc(per_mask, ibc) BIND(C, NAME='BMG_get_bc')

  use ModInterface
  implicit none

  include 'BMG_parameters_f90.h'

  integer(C_INT), value :: per_mask
  integer(C_INT) :: ibc

  integer(C_INT) :: bcmap(0:7)

  bcmap(b'000') = BMG_BCs_definite
  bcmap(b'001') = BMG_BCs_def_per_x
  bcmap(b'010') = BMG_BCs_def_per_y
  bcmap(b'011') = BMG_BCs_def_per_xy
  bcmap(b'100') = BMG_BCs_def_per_z
  bcmap(b'101') = BMG_BCs_def_per_xz
  bcmap(b'110') = BMG_BCs_def_per_yz
  bcmap(b'111') = BMG_BCs_def_per_xyz

  ibc = bcmap(per_mask)

end subroutine BMG_get_bc
