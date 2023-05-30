subroutine BMG_get_bc(per_mask, ibc) BIND(C, NAME='BMG_get_bc')

  use ModInterface
  implicit none

  include 'BMG_parameters_f90.h'

  integer(C_INT), value :: per_mask
  integer(C_INT) :: ibc

  integer(C_INT) :: bcmap(0:7)

  bcmap(0) = BMG_BCs_definite
  bcmap(1) = BMG_BCs_def_per_x
  bcmap(2) = BMG_BCs_def_per_y
  bcmap(3) = BMG_BCs_def_per_xy
  bcmap(4) = BMG_BCs_def_per_z
  bcmap(5) = BMG_BCs_def_per_xz
  bcmap(6) = BMG_BCs_def_per_yz
  bcmap(7) = BMG_BCs_def_per_xyz

  ibc = bcmap(per_mask)

end subroutine BMG_get_bc
