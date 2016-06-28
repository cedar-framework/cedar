module diffusion
  contains

    real(real_t) function dr(i, j, hx, hy)

      use ModInterface
      implicit none

      integer(len_t) :: i, j
      real(real_t) :: hx, hy

      real(real_t) :: xh

      xh = hy/hx

      dr = 1./xh

    end function dr


    real(real_t) function da(i, j, hx, hy)

      use ModInterface
      implicit none

      integer(len_t) :: i, j
      real(real_t) :: hx, hy

      real(real_t) :: yh

      yh = hx/hy

      da = 1./yh

    end function da

end module
