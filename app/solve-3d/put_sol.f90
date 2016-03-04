subroutine sol(i, j, k, hx, hy, hz, val)
  use ModInterface
  implicit none

  integer(len_t) :: i, j, k
  real(real_t) :: hx, hy, hz, x, y, z, val
  real(real_t), parameter :: pi=4*atan(1.d0)

  x = i*hx
  y = j*hy
  z = k*hz

  val = sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)

end subroutine sol


subroutine put_sol(q, ii, jj, kk,&
     hx, hy, hz) BIND(C, NAME='put_sol')
  use ModInterface
  implicit none

  integer(len_t), value :: ii, jj, kk
  real(real_t) :: q(ii,jj,kk)
  real(real_t), value :: hx, hy, hz

  integer :: i, i1, j, j1, k, k1
  real(real_t) :: xh, yh, zh, val

  i1=ii-1
  j1=jj-1
  k1=kk-1

  do k=2,k1
     do j=2,j1
        do i=2,i1
           call sol(i,j,k,hx,hy,hz,val)
           q(i,j,k) = val
        enddo
     enddo
  enddo

end subroutine put_sol
