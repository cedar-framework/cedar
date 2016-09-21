module diffusion
  contains

    real(real_t) function dr(i, j, hx, hy)

      use ModInterface
      implicit none

      integer(len_t) :: i, j
      real(real_t) :: hx, hy

      real(real_t) :: xh

      dr = 1.

    end function dr


    real(real_t) function da(i, j, hx, hy)

      use ModInterface
      implicit none

      integer(len_t) :: i, j
      real(real_t) :: hx, hy

      real(real_t) :: yh

      da = 1.

    end function da

end module
