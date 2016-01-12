subroutine rhs(i, j, k, hx, hy, hz, fs)
  use ModInterface
  implicit none

  integer(len_t) :: i, j, k
  real(real_t) :: hx, hy, hz, x, y, z, fs
  real(real_t), parameter :: pi=4*atan(1.d0)

  x = i*hx
  y = j*hy
  z = k*hz

  fs = 12*(pi**2) * sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)

end subroutine rhs


subroutine putf(so, qf,&
     ii, jj, kk,&
     hx, hy, hz) BIND(C, NAME='putf')

  use ModInterface
  implicit none

  integer(len_t), value :: ii, jj, kk
  real(real_t) :: so(ii,jj,kk,4), qf(ii,jj,kk)
  real(real_t), value :: hx, hy, hz

  integer :: kb, kp, kpw, kps
  parameter ( kp=1, kpw=2, kps=3, kb=4 )

  integer :: i, i1, i2, j, j1, j2, k, k1
  real(real_t) :: xh, yh, zh, h2, fs

  h2=hx*hy*hz
  xh=hy*hz/hx
  yh=hx*hz/hy
  zh=hx*hy/hz
  i1=ii-1
  j1=jj-1
  i2=ii-2
  j2=jj-2
  k1=kk-1

  do k=2,k1
     do j=3,j1
        do i=2,i1
           so(i,j,k,kps) = 1.0*yh;
        enddo
     enddo
  enddo

  do k=2,k1
     do j=2,j1
        do i=3,i1
           so(i,j,k,kpw)= 1.0*xh;
        enddo
     enddo
  enddo

  do k=3,k1
     do j=2,j1
        do i=2,i1
           so(i,j,k,kb) = 1.0*zh;
        enddo
     enddo
  enddo

  do k=2,k1
     do j=2,j1
        do i=2,i1
           so(i,j,k,kp) = so(i,j+1,k,kps) + so(i,j,k+1,kb) + so(i,j,k,kpw)&
                + so(i+1,j,k,kpw) + so(i,j,k,kb) + so(i,j,k,kps)
           call rhs(i,j,k,hx,hy,hz,fs)
           qf(i,j,k) = fs*h2
        enddo
     enddo
  enddo

end subroutine putf
