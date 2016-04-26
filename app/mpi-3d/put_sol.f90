subroutine sol(i, j, k, hx, hy, hz, val)
  use ModInterface
  implicit none

  integer(len_t) :: i, j, k
  real(real_t) :: hx, hy, hz, x, y, z, val
  real(real_t), parameter :: pi=4*atan(1.d0)

  x = (i-1)*hx
  y = (j-1)*hy
  z = (k-1)*hz

  val = sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)

end subroutine sol


subroutine put_sol(q, nlx, nly, nlz,&
     igs, jgs, kgs,&
     hx, hy, hz) BIND(C, NAME='put_sol')
  use ModInterface
  implicit none

  integer(len_t), value :: nlx, nly, nlz
  integer(len_t), value :: iGs, jGs, kGs

  real(real_t) :: q(nlx+2,nly+2,nlz+2)
  real(real_t), value :: hx, hy, hz

  integer(len_t) :: i, i1, j, j1, k, k1
  integer(len_t) :: is, js, kss
  real(real_t) :: xh, yh, zh, val

  i1=nlx+1
  j1=nly+1
  k1=nlz+1

  do k=2,k1
     do j=2,j1
        do i=2,i1
           is = iGs+i-1
           js = jGs+j-1
           kss = kGs+k-1

           call sol(is,js,kss,hx,hy,hz,val)
           q(i,j,k) = val
        enddo
     enddo
  enddo

end subroutine put_sol
